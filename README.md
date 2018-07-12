# M_hotspot
M_hotspot
## author: swyeh

## date: 20180708

## R. version: R version 3.4.2 (2017-09-28)

## used packages: grid, ggplot2, ggrepel, shape, e1071, plotrix, plyr



M_hotspot {swyeh defined}



The core function of M-hotspot

Descritpion
User can find gene mutation hotspot with M-hotspot function.

Usage

M_hotspot (input, total_prop, compare_hotspot, patient_density, hotspot_density, alph1, beta,
          peak_region_min, peak_region_max, remain_gap, num_center, outputpdf,
          file_output, filter_output)

Arguments

input: The table with the somatic mutation amino acid position.
columns definition：(fixed header names)
V1: ENSEMBL ID, GNSP.
V2: Amino acid position of somatic mutation.
V3: TCGA patient's S_ID.
V4: Gene symbol.
V5: Cancer type.

total_prop: Numeric.one of the condition filters. Number of the hotspot mutations > number of gene mutations 5%.
compare_hotspot: Numeric.one of the condition filters. the hotspot with No. mutations > 20% mutations of the strongest hotspot.
patient_density: Numeric. one of the condition filters. The hotspot density per hundred people > 10%.
hotspot_density:Numeric.one of the condition filters. The hotspot density > 2 (each hotspot mutations / each hotspot length)

alph1: Numeric. Parameter of density function. Default is 2.5
beta: Numeric. Parameter of mountain subtractive clustering function. Default is 0.005.

peak_region_min: Numeric. The minimum of the gap criteria. Default is 3.
peak_region_max: Numeric. The maximum of the gap criteria. Default is 10.
remain_gap: Numeric. The top 20% small gap’s information.
num_center: Numeric. Deciding how many hotspot centers to find. Default is 10.

output_pdf: Images output path.
file_output: Hotspots details output path.
filter_output: Hotspots filter details output path.

Value (Output)
	Images, save as pdf file.
  Character, save as txt file.
