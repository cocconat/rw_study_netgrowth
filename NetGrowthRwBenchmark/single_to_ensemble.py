import numpy as np
import matplotlib.pyplot as plt

def SplitSwcFile(input_file, element_type):
    """
    From a single neuron swc file select one element type and cut the tree in a set of
    segments without branching.
    Return them in a list
    """

    f = open(input_file, 'r')
    segments=[[]]
    for line in f:
        if not line.startswith('#') and line.strip():
            if int(line.split()[1]) in element_type:
                sample_number= int(line.split()[0])
                parent_sample= int(line.split()[-1])
                if parent_sample == sample_number-1:
                    segments[-1].append(line)
                else:
                    segments.append([])

    return segments

def SegmentToPath(segment):
    """
    Convert the segment from SWC file into path,
    The enriched information will contain the angle difference between two consecutive pieces, the angle is required to clean the path from
    sudden curves.
    """
    matrix = np.ndarray((3,len(segment)))
    x_0 = float(segment[0].split()[2])
    y_0 = float(segment[0].split()[3])
    x = float(segment[1].split()[2])
    y = float(segment[1].split()[3])
    theta_0 = np.arctan2(-(y-y_0),-(x-x_0))
    for n, line in enumerate(segment):
        x = float(line.split()[2])
        y = float(line.split()[3])
        theta = np.arctan2(-(y-y_0),-(x-x_0))
        matrix[0][n]= x
        matrix[1][n]= y
        matrix[2][n]= theta - theta_0
        x_0 = x
        y_0 = y
        theta_0 = theta
    # plt.plot(matrix[0,:],matrix[1,:])
    # print(matrix[2,matrix[2,:]>1])
    # print(matrix[2,:])
    # plt.scatter(matrix[0,matrix[2,:]>1],matrix[1,matrix[2,:]>1], c='g')
        # if matrix[2][n] > 0.5:
            # return matrix, segment[n:]
    return matrix

def SegmentFromThresh(segment,thresh):
    """
    Cut the segment into subsegments when the angle is bigger
    than threshold
    """
    if np.any(np.abs(segment[2,:])>thresh):
        breakpoints = np.where(np.abs(segment[2,:])>thresh)[0][1:]
        broken=[]
        start=0
        stop=0
        for stop in breakpoints:
            broken.append(segment[:,start:stop-1])
            start=stop
        broken.append(segment[:,stop:])
        return broken
    else:
        return [segment]

def omologate(segments, thresh):
    """
    cut all the segments to 'thresh' length and move each one to start in (0,0)
    remove angle information too
    """
    paths=[]
    for n, segment in enumerate(segments):
        if segment.shape[1]> thresh:
            paths.append((np.subtract(segment.transpose() , segment[:,0]).transpose())[:2,:thresh])
    return paths

def SwcToSegments(input_file, angle_thresh, length_thresh, element_type=[3]):
    segments = SplitSwcFile(input_file, element_type)
    paths=[]
    for seg in segments:
        if len(seg)>2:
            segment = SegmentToPath(seg)
            paths.extend(SegmentFromThresh(segment,angle_thresh))
    paths =np.array(omologate(paths, length_thresh))
    if paths.shape[0] ==0:
        raise ValueError("Segments list is empty")
    return paths

def SegmentsToNetgrowth(paths, name, info):
    """
    NetGrowth_data:
    {"neurons":
        {
        1:{"1":neuronID, "data":swc array}
        2:{"2":neuronID, "data":swc array}

        n:{"3":neuronID, "data":swc array}
        }
    , "info":simulation_parameters}
    """
    NetGrowth_data ={}
    NetGrowth_data["info"]={"name":name,"description":info}
    neurons = {}
    for n, path in enumerate(paths):
        neurons[n]={"gid":n, "data":path.transpose()}
    NetGrowth_data["neurons"]=neurons
    return NetGrowth_data


# for n,x in enumerate(paths):
    # plt.plot(x[0,:],x[1,:])

# plt.show()
