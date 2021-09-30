# Thermophysical Coefficients
**By Keivan A.Gh.**


## Details
As explained in *Chapter 4.8* of the [report](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/Documents/Master's%20Thesis%20-%20Final%20Edition.pdf), the [thermo.compressibleGas](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/OpenFOAM%20v.7/New%20Thermophysical%20Coefficients/thermo.compressibleGas) dictionary contains the last version of the newly generated thermophysical coefficients for hydrogen, oxygen, and nitrogen molecules. These coefficients are required for the implemented [nasa9Poly](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/src/thermophysicalModels/specie/thermo/nasa9Poly) thermodynamic model and [nasaPoly](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/src/thermophysicalModels/specie/transport/nasaPoly) transport model. These coefficients can be used from the boiling points of these species up to 6000 K.

However, in the case of using these two thermophysical models for other temperature intervals or even other species, you can use the original [NASA-9](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/OpenFOAM%20v.7/New%20Thermophysical%20Coefficients/2002_(NASA%209)_McBride%2C%20Gordon_NASA%20Glenn%20Coefficients%20for%20Calculating%20Thermodynamic%20Properties%20of%20Individual%20Species.pdf) coefficients for the *nasa9Poly* thermodynamic model and the original [NASA-7](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/OpenFOAM%20v.7/New%20Thermophysical%20Coefficients/1993_(NASA%207)_McBride%2C%20Gordon_Coefficients%20for%20Calculating%20Thermodynamic%20and%20Transport%20Properties%20of%20Individual%20Species.pdf) coefficients for the *nasaPoly* transport model.  


## Manual
In order to use the [thermo.compressibleGas](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/OpenFOAM%20v.7/New%20Thermophysical%20Coefficients/thermo.compressibleGas) dictionary, you can download the file and move it to the *constant* folder of your OpenFOAM case and add the below code line into your *thermophysicalProperties* file.
```
foamChemistryThermoFile "$FOAM_CASE/constant/thermo.compressibleGas";
```
In the same file, you should also change the thermodynamic and transport models into *nasa9Poly* and *nasaPoly*.

## Citation
For citation, please refer to the [README](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/README.md) file in the main folder.
