This folder contains 1D, 2D and 3D functions for computing sums within an arbitrarily large rectangular ROI in an image from its integral array in constant time. 

In particular, `IntegralSum.jl` contains the most basic functions for computing rectangular-ROI sums from integral arrays. These functions (`IntegralSum(...)` and `IntegralSum_unsafe(...)`) accept four arguments: 
    * the integral array
    * the top-left-front corner of the ROI of interest
    * the bottom-right-corner of the ROI of interest
    * an optional factor, in case you wish to multiply the sums by some value. It defaults to 1.

The functions in `IntegralSumN.jl` combine "integral sums" and "the number of elements in the ROI". In particular, the average value within a ROI can be computed by dividing an integral sum by the area of the ROI. This is offered by `IntegralAvg(...)` and `IntegralAvg_unsafe(...)`. Perhaps less useful is to multiply the integral sum by the area of the ROI, which is offered by `IntegralSumN(...)` and `IntegralSumN_unsafe(...)`. 

The files `IntegralSum_ring.jl` and `IntegralSumN_ring.jl` are a small example of an interesting possiblity with integral arrays, which is to compute sums within "non-rectangular shapes". This can be achieved by decomposing the shape of interest into rectangular regions, and combining integral sums (through addition and substraction) in a somewhat clever way to end up with the sums within the shape of interest. These two files allow to compute the sums within a ring by substracting the sums in small rectangle from the sums in a larger rectangle. 