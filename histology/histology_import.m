%% Import downsampled stains
%{
These are the stains that our collaborator performed. The grayscale images
can be directly used to compute the stain density surrounding the EPVS
because there is no counterstain. The deconvolved images had a
counterstain and therefore had to be deconvolved to extract the target
stain.

Purpose of this script:
- import stains
    - convert to uint8 grayscale (if not already)
- import EPVS annotations
    - convert to logical
- import tissue border
    - convert to logical
- systematic control measurements:
    - evenly-spaced control measurements (exclude EPVS)
- create struct of results

Grayscale Images:
- LHE
- Gallyas

Deconvolved:
- CD68
- Fibrin
- GFAP

%}


%% LHE


%% Gallyas

%% CD68

%% Fibrin

%% GFAP