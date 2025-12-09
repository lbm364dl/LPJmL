# Run LPJmL in a subset of cells

- For adjacent cell range use options `startgrid` and `endgrid` in `lpjml_config.cjson`
- For arbitrary subset of cells there's a workaround. In the [soil types input](#soil-types), if the cell contain a null value, LPJmL will skip that cell. So you can create a soil types input file where only the cells you want to run have a valid soil type, and the rest are null values. [Reference](https://github.com/PIK-LPJmL/LPJmL/issues/71)

# Currently used LPJmL inputs

### Soil types

#### `input.cjson`:

```json
  "soil" :         { "id" : 41, "fmt" : "cdf", "var": "soil_type", "name" : "soil/soil_30arcmin_13_types.nc"},
```

#### Layers

Static, no time series, 1 layer.

#### Values

Types of soil defined in `lpjml_config.cjson`:

```json
  "soilmap" : [null,"clay", "silty clay", "sandy clay", "clay loam", "silty clay loam",
              "sandy clay loam", "loam", "silt loam", "sandy loam", "silt",
              "loamy sand", "sand", "rock and ice"],
```

Raster values are numbers from 1 to 13, indicating the index in the above list (starting from 0). Careful with using the number 0. If you want to skip a cell remember to use `null` directly. In terra R package, you can use `NA` to represent null values.

### Soil pH 

#### `input.cjson`:

```json
  "soilpH" :        { "id" : 26, "fmt" : "clm",  "name" : "soil/soil_pH_30arcmin.bin"},
```

#### Layers

Static, no time series, 1 layer.

#### Values

Decimal values between 0 and 14 representing the pH of the soil in each grid cell.


### Land fraction

#### `input.cjson`:

```json
  "landfrac" :     { "id" : 44, "fmt" : "clm",  "name" : "gadm/landfrac_gadm_30arcmin.bin"},
```

#### Layers

Static, no time series, 1 layer.

#### Values

Decimal values between 0 and 1 representing the fraction of land in each grid cell. Value 0 means no land, value 1 means full land. You can have intermediate values like 0.3 meaning 30% of the grid cell is land.

### Lakes and rivers fraction

#### `input.cjson`:

```json
  "lakes" :     { "id" : 36, "fmt" : "clm",  "name" : "lakes_rivers/glwd_lakes_and_rivers_30arcmin.bin"},
```

#### Layers

Static, no time series, 1 layer.

#### Values

Decimal values between 0 and 1 representing the fraction of water bodies in each grid cell. Here water bodies include lakes and rivers. Value 0 means no water bodies, value 1 means full water bodies. You can have intermediate values like 0.2 meaning 20% of the grid cell is covered by lakes and rivers.

### Near-surface temperature

#### `input.cjson`:

```json
  "temp" : { "fmt" : "cdf", "var" : "tmp", "name" : "climate/cru_ts_3_10.1901.2009.tmp.dat.nc"},
```

#### Layers

Monthly/daily time series. One layer per month/day. Set the dates properly for each layer, e.g., named "1901-01-01", "1901-02-01", etc. for monthly data. The day of the month is not important for monthly data, just the year and month, but the date must be in that format.

#### Values

Decimal values in degrees Celsius.


### Precipitation

#### `input.cjson`:

```json
  "prec" : { "fmt" : "cdf", "var" : "pre", "name" : "climate/cru_ts_3_10_01.1901.2009.pre.dat.nc"},
```

#### Layers

Monthly/daily time series. One layer per month/day. Set the dates properly for each layer, e.g., named "1901-01-01", "1901-02-01", etc. for monthly data. The day of the month is not important for monthly data, just the year and month, but the date must be in that format.

#### Values

Decimal values in mm (per day or month depending on the layers).


### Cloud cover

#### `input.cjson`:

```json
  "cloud" : { "fmt" : "cdf", "var" : "cld", "name" : "climate/cru_ts_3_10.1901.2009.cld.dat.nc"},
```

#### Layers

Monthly/daily time series. One layer per month/day. Set the dates properly for each layer, e.g., named "1901-01-01", "1901-02-01", etc. for monthly data. The day of the month is not important for monthly data, just the year and month, but the date must be in that format.

The dataset we downloaded from ISIMIP came with 4 x time layers, where one layer contains the main cloud cover data and the others are metadata about the estimation. Having a single variable called "cld" should be enough, LPJmL will ignore the rest.

#### Values

Decimal values between 0 and 1 representing the fraction of cloud cover. The value 0 means no clouds, the value 1 means full cloud cover.

### Wind speed

#### `input.cjson`:

```json
  "wind" : { "fmt" : "cdf", "var" : "wind", "name" : "climate/wind_gswp3-w5e5_1901_2016_monthly.nc"},
```

#### Layers

Monthly/daily time series. One layer per month/day. Set the dates properly for each layer, e.g., named "1901-01-01", "1901-02-01", etc. for monthly data. The day of the month is not important for monthly data, just the year and month, but the date must be in that format.

#### Values

Decimal values in meters per second

### Wet days

#### `input.cjson`:

```json
  "wetdays" :      { "id" : 12, "fmt" : "clm",  "name" : "CRUDATA_TS3_23/gpcc_v7_cruts3_23_wet_1901_2013.clm"},
```

#### Layers

Monthly time series. One layer per month. Set the dates properly for each layer, e.g., named "1901-01-01", "1901-02-01", etc. for monthly data. The day of the month is not important for monthly data, just the year and month, but the date must be in that format.

#### Values

Decimal values representing the average number of wet days in that month, so should be between 0 and 31.

### NO3- deposition

#### `input.cjson`:

```json
  "no3deposition" : { "id" : 17, "fmt" : "cdf", "var": "noy", "name" : "nitrogen/ndep_noy_histsoc_monthly_1860_2016.nc4"},
```

#### Layers

Monthly time series. One layer per month. Set the dates properly for each layer, e.g., named "1901-01-01", "1901-02-01", etc. for monthly data. The day of the month is not important for monthly data, just the year and month, but the date must be in that format.

#### Values

Decimal values in grams of nitrogen per square meter per month.


### NH4+ deposition

#### `input.cjson`:

```json
  "nh4deposition" : { "id" : 16, "fmt" : "cdf", "var": "nhx", "name" : "nitrogen/ndep_nhx_histsoc_monthly_1860_2016.nc4"},
```

#### Layers

Monthly time series. One layer per month. Set the dates properly for each layer, e.g., named "1901-01-01", "1901-02-01", etc. for monthly data. The day of the month is not important for monthly data, just the year and month, but the date must be in that format.

#### Values

Decimal values in grams of nitrogen per square meter per month.

### CO2 concentration

#### `input.cjson`:

```json
  "co2" :          { "id" : 5, "fmt" : "txt",  "name" : "climate/historical_CO2_annual_1765_2018.txt"},
```

#### Layers

Not a raster, just a time series with one value per year in a tab separated file. Example:

```txt
1790	281.60
1791	281.75
1792	281.89
1793	282.03
1794	282.17
1795	282.30
1796	282.43
1797	282.55
1798	282.67
1799	282.79
1800	282.90
```

### Land use types

#### `input.cjson`:

```json
  "landuse" :      { "id" : 6, "fmt" : "clm",  "name" : "landuse/cft_default_cft_aggregation_30min_1500-2017.bin"},
```

#### Layers

Static, no time series, 1 layer.

#### Values

Types of soil defined in `lpjml_config.cjson`:

```json
  "soilmap" : [null,"clay", "silty clay", "sandy clay", "clay loam", "silty clay loam",
              "sandy clay loam", "loam", "silt loam", "sandy loam", "silt",
              "loamy sand", "sand", "rock and ice"],
```

Raster values are numbers from 1 to 13, indicating the index in the above list (starting from 0). Careful with using the number 0. If you want to skip a cell remember to use `null` directly. In terra R package, you can use `NA` to represent null values.

#### Values

Decimal values representing the CO2 concentration in parts per million (ppm).

### Countries and subdivisions

#### `input.cjson`:

```json
  "countrycode" :  { "id" : 25, "fmt" : "clm",  "name" : "gadm/cow_gadm_30arcmin.bin"},
```

#### Layers

Static, no time series, 2 layers.

Layer 1: country code of the one with the largest area in that grid cell.
Layer 2: subregion code of the one with the largest area in that grid cell.

#### Values

Integer codes representing countries and subregions. 

- The country codes can be found as defined in `include/managepar.h`. 
- TODO: find about the subregion codes. Currently probably not used?

### Land-use fractions

#### `input.cjson`:

```json
  "landuse" :      { "id" : 6, "fmt" : "cdf", "var" : "landuse",  "name" : "landuse/cft_default_cft_aggregation_30min_1900-2017.nc"},
```

#### Layers

Yearly time series. 32 layers (16 rainfed, 16 irrigated). Each 16 layers come in the order defined by the map:

```json
  /* the following array describes the order of the bands in the land use file */
  "landusemap" : ["temperate cereals","rice", "maize", "tropical cereals", "pulses",
                  "temperate roots", "tropical roots", "oil crops sunflower",
                  "oil crops soybean", "oil crops groundnut", "oil crops rapeseed",
                  "sugarcane","others","grassland","biomass grass","biomass tree"],
```

TODO: this is currently not working in NetCDF format. Working on it https://github.com/PIK-LPJmL/LPJmL/issues/73

#### Values

Decimal values between 0 and 1 representing percentage of land use in that grid cell for that crop type and irrigation system.


# Other things about inputs

### Number of bands x2 or x4 for irrigation

Whenever asked for doubled number of bands, it should mean the first half represents rainfed and the second half irrigated. If asked for quadruple number of bands, it should be in this order: no irrigation, surface irrigation, sprinkler, drip.

### Land use input file warning

`WARNING024: Land-use input does not include irrigation systems, suboptimal country values are used.`

Just a warning, so the simulation works anyway. Just keep in mind for possible fix.

Land use input file `cft_default_cft_aggregation_30min_1500-2017` was obtained from LandInG repository. Check if something wrong there. Otherwise it might just be calculated like this, without quadruple bands for irrigation systems.

Suboptimal country values refer to the values of `default_irrig_system` for each country in `par/manage_irrig_systems_with_dummy_laimax_data.cjson`. The four irrigation systems are: `NOIRRIG, SURF, SPRINK, DRIP`.


### Some land use inputs and how to skip them

#### residues_madrat_1850-2015_16bands

16 bands that come from a predefined list in lpjml_config.cjson

```json
  "fertilizermap" : ["temperate cereals","rice", "maize", "tropical cereals", "pulses",
                    "temperate roots", "tropical roots", "oil crops sunflower",
                    "oil crops soybean", "oil crops groundnut", "oil crops rapeseed",
                    "sugarcane","others","grassland","biomass grass","biomass tree"],
```

Value: ratio (between 0 and 1) that represents, from total generated residue, which percentage is kept in the land

Use: https://github.com/PIK-LPJmL/LPJmL/blob/383a0c510153a0a9f1b7a29f2694cfc1476a919d/src/crop/harvest_crop.c#L51

Skippable: look for this line in lpjml_config.cjson and change to "no_residue_remove"

  "residue_treatment" : "read_residue_data", /* residue options: "read_residue_data" (requires input dataset on fraction to remain on field), "no_residue_remove" (all non-harvested plant components to litter), "fixed_residue_remove" (fraction of non-harvested above ground biomass to remain in soil according to param residues_in_soil) */

For the other option of fixed_residue_remove see value e.g. in https://github.com/PIK-LPJmL/LPJmL/blob/383a0c510153a0a9f1b7a29f2694cfc1476a919d/par/lpjparam.cjson#L95

#### tillage_1973_1974-2010_1973_as_till_everywhere.clm

Single band.

Value: 0 for no tillage, 1 for tillage

Skippable: look for this line in lpjml_config.cjson and change to "no" / "all"

```json
  "tillage_type" : "read",              /* Options: "all" (all agr. cells tilled), "no" (no cells tilled) and "read" (tillage dataset used) */
```

#### grassland_lsuha_gswp3_2000-2000.clm

Skippable: look for this line in lpjml_config.cjson and change to false

```json
  "prescribe_lsuha" : true,             /* prescribe livestock unit per hectare, input required (needed for "grazing" = "livestock") */
```

#### sdate_24_bands.clm

Skippable: look for this line in lpjml_config.cjson and change

```json
  "sowing_date_option" : "prescribed_sdate",   /* options: no_fixed_sdate, fixed_sdate, prescribed_sdate */
```

Using "no_fixed_sdate" lets the model dynamically find the sowing dates each year.
Using "fixed_sdate" uses calculated sowing dates from a fixed year, for all years. It comes from this line in lpjml_config.cjson:

```json
  "sdate_fixyear" : 1970,               /* year in which sowing dates shall be fixed, when using fixed_sdate */
```

#### phusum_gswp3-w5e5_obsclim_24_bands_annual_1901-2019.clm

Skippable: look for this line in lpjml_config.cjson and change

```json
  "crop_phu_option" : "prescribed",     /* defines setting for phenological heat unit and vernalization requirements for crops; options: "old", "new", "prescribed" */
```

Not sure what "old" and "new" do but probably should use "new" if it's an improved calculation...

Use: maybe check https://github.com/PIK-LPJmL/LPJmL/blob/master/src/crop/phenology_crop.c

#### fert_N_default_cft_aggregation_30min_1860-2017.clm

Skippable: look for this line in lpjml_config.cjson and change to true

```json
  "fix_fertilization" : false,          /* fertilizer and manure input setting, options: "false" (prescribed time/cell/cft specific input), "true" (global constant rates from lpjmlparam.cjson) */
```
