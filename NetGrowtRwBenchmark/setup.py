import os

from setuptools import setup


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as buf:
        return buf.read()


conf = dict(
        name='rw_corr',
        version='0.1',
        author='cocco',
        license='AGPL',
        packages=['rw_corr'],

        zip_safe=False,
        entry_points={'console_scripts': [
            'rw_corr=rw_corr.main:main',
        ]},
        classifiers=[
          "License :: OSI Approved :: GNU Affero General Public License v3",
          "Operating System :: POSIX :: Linux",
          "Programming Language :: Python :: 3",
        ])


if __name__ == '__main__':
    setup(**conf)
