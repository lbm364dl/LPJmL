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

### [Input name]
#### `input.cjson`:

```json
  "[input_key]" : { "id" : [id], "fmt" : "[format]", "var": "[variable]", "name" : "[path/to/file]" },
```

#### Layers

[Static or temporal, number of layers/bands, time series info.]

#### Values

[Description of value types, units, valid ranges, and any mapping or typology.]



# Other things about inputs

### Number of bands x2 or x4 for irrigation

Whenever asked for doubled number of bands, it should mean the first half represents rainfed and the second half irrigated. If asked for quadruple number of bands, it should be in this order: no irrigation, surface irrigation, sprinkler, drip.

### Landuse input file warning

`WARNING024: Land-use input does not include irrigation systems, suboptimal country values are used.`

Just a warning, so the simulation works anyway. Just keep in mind for possible fix.

Landuse input file `cft_default_cft_aggregation_30min_1500-2017` was obtained from LandInG repository. Check if something wrong there. Otherwise it might just be calculated like this, without quadruple bands for irrigation systems.

Suboptimal country values refer to the values of `default_irrig_system` for each country in `par/manage_irrig_systems_with_dummy_laimax_data.cjson`. The four irrigation systems are: `NOIRRIG, SURF, SPRINK, DRIP`.


### Some landuse inputs and how to skip them

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
