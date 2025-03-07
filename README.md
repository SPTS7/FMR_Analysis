# Magnetic Field Fitting and Analysis

## Overview

This script performs data analysis for magnetic field experiments by fitting a combination of Lorentzian models to the experimental data. It calculates key material properties such as peak frequency, full width at half maximum (FWHM), damping values, and gyromagnetic ratio using different anisotropy models. The script generates plots to visualize the results and exports the fitting results and derived parameters into CSV files for further analysis.

## Dependencies

- `scipy`
- `numpy`
- `matplotlib`
- `pyqtgraph`
- `lmfit`
- `tqdm`

These dependencies can be installed using `pip`:

```bash
pip install scipy numpy matplotlib pyqtgraph lmfit tqdm
```

## File Format

The script expects the input data files to be in `.txt` format. Each file should contain two columns of data:

- **Column 1**: Magnetic field (Oe)
- **Column 2**: Measured signal (arbitrary units)

The files should be named in such a way that the magnetic field value can be extracted from the filename (e.g., `file_100.txt` where `100` is the field value).

## Usage

1. **Place your `.txt` data files** in the same directory as this script.

2. **Configure the script settings**:

   - Set the desired **anisotropy type** by modifying the `Anisotropy` variable. Options are:
     - `"IP"`: In-Plane anisotropy
     - `"OP"`: Out-of-Plane anisotropy
     - `"U110"`: Ultrathin 110 anisotropy
   - Adjust material parameters (e.g., `G`, `Ms`, `K2`, `alpha`, etc.) according to your specific setup.

3. **Run the script**. The script will:

   - Automatically find all `.txt` files in the directory.
   - Extract the magnetic field values from the filenames.
   - Fit Lorentzian models to the data to find peaks and calculate material properties.
   - Generate plots and save them to PNG files.
   - Export the fitting results and derived parameters to CSV files.

4. **Output**:
   - **Plots**: The script generates a series of plots showing the fitted peak, FWHM, damping, and other properties. These are saved as PNG files with the naming pattern `Fit_Peak{peak_number}.png`.
   - **CSV Files**:
     - `fields_fwhm_peak_{peak_number}.csv`: Contains the magnetic field, peak frequency, and FWHM values.
     - `Dampings{peak_number}.csv`: Contains the damping values and related data.
     - `Fit_Values_Peak{peak_number}.csv`: Contains the detailed fit report from `Kittel` and damping models.

## Example Output

1. **Generated Plots**:

   - Magnetic Field vs Peak Frequency
   - Magnetic Field vs FWHM
   - Peak Frequency vs dH
   - Damping vs Magnetic Field

2. **CSV Output** (sample data structure):

```csv
# fields_fwhm_peak_1.csv
Field, Peak, FWHM
100, 2.5, 0.1
200, 2.7, 0.12
...

# Dampings1.csv
x1_Peak, y1_dH, x2_newpeak, y2_dH_Fit, x3_Field, y3_dHa, y3_baseline
2.5, 0.15, 2.52, 0.16, 100, 0.18, 0.2
...

# Fit_Values_Peak1.csv
Kittel Fit Report:
[Detailed Kittel Fit Report]
Damping Fit Report:
[Detailed Damping Fit Report]
```

## Key Functions

- `Lorexponential(x, simetricL, asymetricL, center, sigma, C)`: Defines a combined symmetric and asymmetric Lorentzian model.
- `make_model(num)`: Creates the model for each detected peak.
- `KittelIP(field, G, Hintrinsic)`: Defines the Kittel fitting function for In-Plane anisotropy.
- `KittelOP(field, G, Hintrinsic)`: Defines the Kittel fitting function for Out-of-Plane anisotropy.
- `Ultrathin110(field, G, K2, K4, Ms)`: Defines the Kittel fitting function for Ultrathin 110 anisotropy.
- `damping(x, dH0, G, alpha)`: Defines the damping model.

## Contributing

Feel free to fork this repository and contribute by opening issues or submitting pull requests for improvements. Make sure to follow the existing coding style and conventions.

## License

This project is licensed under the MIT License.
