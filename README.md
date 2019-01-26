# wavelet-networks

This project is an implementation of the paper of Qinghua Zhang and Albert Benveniste on wavelet networks: [link](https://ieeexplore.ieee.org/document/165591/).

## Getting started

```{shell}
> mkdir build
> cd build
> cmake ..
> make
> ./waveletnet
```

In the current state, the `main.cpp` file instantiates a wavelet network trained on noisy sinus data and evaluate the results using smooth data.
The networks contains 11 wavelons and is trained with a learning rate of 0.05.
