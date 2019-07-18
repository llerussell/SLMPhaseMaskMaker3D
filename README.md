# SLMPhaseMaskMaker3D
Generates phase masks for display on an SLM using a weighted GS algorithm (Martin Persson's [HOTlab](https://github.com/MartinPersson/HOTlab))

To calibrate the generation of SLM phase masks that target ROIs in an imaging stack see [SLMTransformMaker](https://github.com/llerussell/SLMTransformMaker3D)

SLM phase mask generation is incorporated into our all-optical experiment control software [Naparm](https://github.com/llerussell/Naparm)

![Imgur](http://i.imgur.com/fNY3ewP.jpg)

# Requirements
* CUDA 7/8

## Example standalone usage
1. Simply run the script and load a targets image (black 512x512 with white pixels corresponding to target centroids)
2. In code:
``` matlab
% Make and save phase mask(s)
[PhaseMask, TransformedSLMTargets] = SLMPhaseMaskMakerCUDA3D(...
    'Points', points,...
    'Save', true,...
    'SaveName', ['ExamplePhaseMask.tif'],...
    'Do2DTransform', false,...
    'Do3DTransform', false,...
    'AutoAdjustWeights', false);

    % where points is a n * 4 array (num of points * [x y z intensity])
```
