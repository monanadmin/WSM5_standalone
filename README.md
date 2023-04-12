# WSM5_standalone
## Single Moment 5-Class (WSM5) Cloud Microphysics - Stand alone package

### Description
The standalone package was created is bases on paper *Improved GPU/CUDA Based Parallel Weather and Research Forecast (WRF) Single Moment 5-Class (WSM5) Cloud Microphysics* [1].

### Abstract

The Weather Research and Forecasting (WRF) model is an atmospheric simulation system which is designed for both operational and research use. WRF is currently in operational use at the National Oceanic and Atmospheric Administration (NOAA)'s national weather service as well as at the air force weather agency and meteorological services worldwide. Getting weather predictions in time using latest advances in atmospheric sciences is a challenge even on the fastest super computers. Timely weather predictions are particularly useful for severe weather events when lives and property are at risk. Microphysics is a crucial but computationally intensive part of WRF. WRF Single Moment 5-class (WSM5) microphysics scheme represents fallout of various types of precipitation, condensation and thermodynamics effects of latent heat release. Therefore, to expedite the computation process, Graphics Processing Units (GPUs) appear an attractive alternative to traditional CPU architectures. In this paper, we accelerate the WSM5 microphysics scheme on GPUs and obtain a considerable speedup thereby significantly reducing the processing time. Such high performance and computationally efήcient GPUs allow us to use higher resolution WRF forecasts. The use of high resolution WRF enables us to compute microphysical processes for increasingly small clouds and water droplets. To implement WSM5 scheme on GPUs, the WRF code was rewritten into CUDA C, a high level data-parallel programming language used on NVIDIA GPU. We observed a reduction in processing time from 16928 ms on CPU to 43.5 ms on a Graphics Processing Unit (GPU). We obtained a speedup of 389× without I/O using a single GPU. Taking I/O transfer times into account, the speedup obtained is 206×. The speedup was further increased by using four GPUs, speedup being 1556× and 357× for without I/O and with I/O, respectively.

### Compiling and running the code

```
$ git clone https://github.com/monanadmin/WSM5_standalone.git
$ cd WSM5_standalone
$ git checkout develop
$ cd datain
$ wget http://ftp.cptec.inpe.br/pesquisa/dmdcc/monan/standalone_codes/WSM5/datain/WSM05_dataIn-03600-00001.bin
$ cd ../scripts
$ ./WSM5_standalone.sh
```

### WSM5 Cloud Microphysics Schema

![Flowchart-of-the-microphysics-processes-in-the-WSM5-scheme](https://user-images.githubusercontent.com/6113640/231156459-dd621b35-0a64-43ac-89c4-6f40b60e7c7a.png)

![image](https://user-images.githubusercontent.com/6113640/231156726-b9a39fd6-24fe-43e1-964d-9e94f26d168c.png)


### References

#### 1. [J. Mielikainen, B. Huang, H. -L. A. Huang and M. D. Goldberg, "Improved GPU/CUDA Based Parallel Weather and Research Forecast (WRF) Single Moment 5-Class (WSM5) Cloud Microphysics," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 5, no. 4, pp. 1256-1265, Aug. 2012, doi: 10.1109/JSTARS.2012.2188780.](https://ieeexplore.ieee.org/document/6182591)
