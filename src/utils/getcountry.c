/**************************************************************************************/
/**                                                                                \n**/
/**                  g  e  t  c  o  u  n  t  r  y  .  c                            \n**/
/**                                                                                \n**/
/**     Extract country from grid file                                             \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#undef USE_MPI
#include "lpj.h"

#define USAGE "Usage: %s [-list] [-longheader] countryfile gridfile outfile country ...\n"

#define NCOUNTRY 257 /* number of countries defined in LPJmL */

typedef struct
{
  int code;     /* LPJmL country code */
  const char *name;
  char abbrev[4]; /* ISO 3166-1 alpha-3 abbreviation for country */
} Countryname;

static Countryname countrynames[NCOUNTRY]=
{
  {Aruba, "Aruba", "ABW"},
  {Afghanistan, "Afghanistan", "AFG"},
  {Angola, "Angola", "AGO"},
  {Anguilla, "Anguilla", "AIA"},
  {Aland, "Åland", "ALA"},
  {Albania, "Albania", "ALB"},
  {Andorra, "Andorra", "AND"},
  {United_Arab_Emirates, "United Arab Emirates", "ARE"},
  {Argentina, "Argentina", "ARG"},
  {Armenia, "Armenia", "ARM"},
  {American_Samoa, "American Samoa", "ASM"},
  {Antarctica, "Antarctica", "ATA"},
  {French_Southern_Territories, "French Southern Territories", "ATF"},
  {Antigua_and_Barbuda, "Antigua and Barbuda", "ATG"},
  {Austria, "Austria", "AUT"},
  {Azerbaijan, "Azerbaijan", "AZE"},
  {Burundi, "Burundi", "BDI"},
  {Belgium, "Belgium", "BEL"},
  {Benin, "Benin", "BEN"},
  {Bonaire_Sint_Eustatius_and_Saba, "Bonaire, Sint Eustatius and Saba", "BES"},
  {Burkina_Faso, "Burkina Faso", "BFA"},
  {Bangladesh, "Bangladesh", "BGD"},
  {Bulgaria, "Bulgaria", "BGR"},
  {Bahrain, "Bahrain", "BHR"},
  {Bahamas, "Bahamas", "BHS"},
  {Bosnia_and_Herzegovina, "Bosnia and Herzegovina", "BIH"},
  {Saint_Barthelemy, "Saint-Barthélemy", "BLM"},
  {Belarus, "Belarus", "BLR"},
  {Belize, "Belize", "BLZ"},
  {Bermuda, "Bermuda", "BMU"},
  {Bolivia, "Bolivia", "BOL"},
  {Barbados, "Barbados", "BRB"},
  {Brunei, "Brunei", "BRN"},
  {Bhutan, "Bhutan", "BTN"},
  {Bouvet_Island, "Bouvet Island", "BVT"},
  {Botswana, "Botswana", "BWA"},
  {Central_African_Republic, "Central African Republic", "CAF"},
  {Cocos_Islands, "Cocos Islands", "CCK"},
  {Switzerland, "Switzerland", "CHE"},
  {Chile, "Chile", "CHL"},
  {Cote_Ivoire, "Côte d'Ivoire", "CIV"},
  {Cameroon, "Cameroon", "CMR"},
  {Democratic_Republic_of_the_Congo, "Democratic Republic of the Congo", "COD"},
  {Republic_of_Congo, "Republic of Congo", "COG"},
  {Cook_Islands, "Cook Islands", "COK"},
  {Colombia, "Colombia", "COL"},
  {Comoros, "Comoros", "COM"},
  {Cape_Verde, "Cape Verde", "CPV"},
  {Costa_Rica, "Costa Rica", "CRI"},
  {Cuba, "Cuba", "CUB"},
  {Curacao, "Curaçao", "CUW"},
  {Christmas_Island, "Christmas Island", "CXR"},
  {Cayman_Islands, "Cayman Islands", "CYM"},
  {Cyprus, "Cyprus", "CYP"},
  {Czech_Republic, "Czech Republic", "CZE"},
  {Germany, "Germany", "DEU"},
  {Djibouti, "Djibouti", "DJI"},
  {Dominica, "Dominica", "DMA"},
  {Denmark, "Denmark", "DNK"},
  {Dominican_Republic, "Dominican Republic", "DOM"},
  {Algeria, "Algeria", "DZA"},
  {Ecuador, "Ecuador", "ECU"},
  {Egypt, "Egypt", "EGY"},
  {Eritrea, "Eritrea", "ERI"},
  {Western_Sahara, "Western Sahara", "ESH"},
  {Spain, "Spain", "ESP"},
  {Estonia, "Estonia", "EST"},
  {Ethiopia, "Ethiopia", "ETH"},
  {Finland, "Finland", "FIN"},
  {Fiji, "Fiji", "FJI"},
  {Falkland_Islands, "Falkland Islands", "FLK"},
  {France, "France", "FRA"},
  {Faroe_Islands, "Faroe Islands", "FRO"},
  {Micronesia, "Micronesia", "FSM"},
  {Gabon, "Gabon", "GAB"},
  {United_Kingdom, "United Kingdom", "GBR"},
  {Georgia, "Georgia", "GEO"},
  {Guernsey, "Guernsey", "GGY"},
  {Ghana, "Ghana", "GHA"},
  {Gibraltar, "Gibraltar", "GIB"},
  {Guinea, "Guinea", "GIN"},
  {Guadeloupe, "Guadeloupe", "GLP"},
  {Gambia, "Gambia", "GMB"},
  {Guinea_Bissau, "Guinea-Bissau", "GNB"},
  {Equatorial_Guinea, "Equatorial Guinea", "GNQ"},
  {Greece, "Greece", "GRC"},
  {Grenada, "Grenada", "GRD"},
  {Greenland, "Greenland", "GRL"},
  {Guatemala, "Guatemala", "GTM"},
  {French_Guiana, "French Guiana", "GUF"},
  {Guam, "Guam", "GUM"},
  {Guyana, "Guyana", "GUY"},
  {Hong_Kong, "Hong Kong", "HKG"},
  {Heard_Island_and_McDonald_Islands, "Heard Island and McDonald Islands", "HMD"},
  {Honduras, "Honduras", "HND"},
  {Croatia, "Croatia", "HRV"},
  {Haiti, "Haiti", "HTI"},
  {Hungary, "Hungary", "HUN"},
  {Indonesia, "Indonesia", "IDN"},
  {Isle_of_Man, "Isle of Man", "IMN"},
  {British_Indian_Ocean_Territory, "British Indian Ocean Territory", "IOT"},
  {Ireland, "Ireland", "IRL"},
  {Iran, "Iran", "IRN"},
  {Iraq, "Iraq", "IRQ"},
  {Iceland, "Iceland", "ISL"},
  {Israel, "Israel", "ISR"},
  {Italy, "Italy", "ITA"},
  {Jamaica, "Jamaica", "JAM"},
  {Jersey, "Jersey", "JEY"},
  {Jordan, "Jordan", "JOR"},
  {Japan, "Japan", "JPN"},
  {Kazakhstan, "Kazakhstan", "KAZ"},
  {Kenya, "Kenya", "KEN"},
  {Kyrgyzstan, "Kyrgyzstan", "KGZ"},
  {Cambodia, "Cambodia", "KHM"},
  {Kiribati, "Kiribati", "KIR"},
  {Saint_Kitts_and_Nevis, "Saint Kitts and Nevis", "KNA"},
  {South_Korea, "South Korea", "KOR"},
  {Kuwait, "Kuwait", "KWT"},
  {Laos, "Laos", "LAO"},
  {Lebanon, "Lebanon", "LBN"},
  {Liberia, "Liberia", "LBR"},
  {Libya, "Libya", "LBY"},
  {Saint_Lucia, "Saint Lucia", "LCA"},
  {Liechtenstein, "Liechtenstein", "LIE"},
  {Sri_Lanka, "Sri Lanka", "LKA"},
  {Lesotho, "Lesotho", "LSO"},
  {Lithuania, "Lithuania", "LTU"},
  {Luxembourg, "Luxembourg", "LUX"},
  {Latvia, "Latvia", "LVA"},
  {Macao, "Macao", "MAC"},
  {Saint_Martin, "Saint-Martin", "MAF"},
  {Morocco, "Morocco", "MAR"},
  {Monaco, "Monaco", "MCO"},
  {Moldova, "Moldova", "MDA"},
  {Madagascar, "Madagascar", "MDG"},
  {Maldives, "Maldives", "MDV"},
  {Mexico, "Mexico", "MEX"},
  {Marshall_Islands, "Marshall Islands", "MHL"},
  {Macedonia, "Macedonia", "MKD"},
  {Mali, "Mali", "MLI"},
  {Malta, "Malta", "MLT"},
  {Myanmar, "Myanmar", "MMR"},
  {Montenegro, "Montenegro", "MNE"},
  {Mongolia, "Mongolia", "MNG"},
  {Northern_Mariana_Islands, "Northern Mariana Islands", "MNP"},
  {Mozambique, "Mozambique", "MOZ"},
  {Mauritania, "Mauritania", "MRT"},
  {Montserrat, "Montserrat", "MSR"},
  {Martinique, "Martinique", "MTQ"},
  {Mauritius, "Mauritius", "MUS"},
  {Malawi, "Malawi", "MWI"},
  {Malaysia, "Malaysia", "MYS"},
  {Mayotte, "Mayotte", "MYT"},
  {Namibia, "Namibia", "NAM"},
  {New_Caledonia, "New Caledonia", "NCL"},
  {Niger, "Niger", "NER"},
  {Norfolk_Island, "Norfolk Island", "NFK"},
  {Nigeria, "Nigeria", "NGA"},
  {Nicaragua, "Nicaragua", "NIC"},
  {Niue, "Niue", "NIU"},
  {Netherlands, "Netherlands", "NLD"},
  {Norway, "Norway", "NOR"},
  {Nepal, "Nepal", "NPL"},
  {Nauru, "Nauru", "NRU"},
  {New_Zealand, "New Zealand", "NZL"},
  {Oman, "Oman", "OMN"},
  {Pakistan, "Pakistan", "PAK"},
  {Panama, "Panama", "PAN"},
  {Pitcairn_Islands, "Pitcairn Islands", "PCN"},
  {Peru, "Peru", "PER"},
  {Philippines, "Philippines", "PHL"},
  {Palau, "Palau", "PLW"},
  {Papua_New_Guinea, "Papua New Guinea", "PNG"},
  {Poland, "Poland", "POL"},
  {Puerto_Rico, "Puerto Rico", "PRI"},
  {North_Korea, "North Korea", "PRK"},
  {Portugal, "Portugal", "PRT"},
  {Paraguay, "Paraguay", "PRY"},
  {Palestina, "Palestina", "PSE"},
  {French_Polynesia, "French Polynesia", "PYF"},
  {Qatar, "Qatar", "QAT"},
  {Reunion, "Reunion", "REU"},
  {Romania, "Romania", "ROU"},
  {Rwanda, "Rwanda", "RWA"},
  {Saudi_Arabia, "Saudi Arabia", "SAU"},
  {Sudan, "Sudan", "SDN"},
  {Senegal, "Senegal", "SEN"},
  {Singapore, "Singapore", "SGP"},
  {South_Georgia_and_the_South_Sandwich_Islands, "South Georgia and the South Sandwich Islands", "SGS"},
  {Saint_Helena, "Saint Helena", "SHN"},
  {Svalbard, "Svalbard and Jan Mayen", "SJM"},
  {Solomon_Islands, "Solomon Islands", "SLB"},
  {Sierra_Leone, "Sierra Leone", "SLE"},
  {El_Salvador, "El Salvador", "SLV"},
  {San_Marino, "San Marino", "SMR"},
  {Somalia, "Somalia", "SOM"},
  {Saint_Pierre_and_Miquelon, "Saint Pierre and Miquelon", "SPM"},
  {Serbia, "Serbia", "SRB"},
  {South_Sudan, "South Sudan", "SSD"},
  {Sao_Tome_and_Principe, "São Tomé and Príncipe", "STP"},
  {Suriname, "Suriname", "SUR"},
  {Slovakia, "Slovakia", "SVK"},
  {Slovenia, "Slovenia", "SVN"},
  {Sweden, "Sweden", "SWE"},
  {Swaziland, "Swaziland", "SWZ"},
  {Sint_Maarten, "Sint Maarten", "SXM"},
  {Seychelles, "Seychelles", "SYC"},
  {Syria, "Syria", "SYR"},
  {Turks_and_Caicos_Islands, "Turks and Caicos Islands", "TCA"},
  {Chad, "Chad", "TCD"},
  {Togo, "Togo", "TGO"},
  {Thailand, "Thailand", "THA"},
  {Tajikistan, "Tajikistan", "TJK"},
  {Tokelau, "Tokelau", "TKL"},
  {Turkmenistan, "Turkmenistan", "TKM"},
  {Timor_Leste, "Timor-Leste", "TLS"},
  {Tonga, "Tonga", "TON"},
  {Trinidad_and_Tobago, "Trinidad and Tobago", "TTO"},
  {Tunisia, "Tunisia", "TUN"},
  {Turkey, "Turkey", "TUR"},
  {Tuvalu, "Tuvalu", "TUV"},
  {Taiwan, "Taiwan", "TWN"},
  {Tanzania, "Tanzania", "TZA"},
  {Uganda, "Uganda", "UGA"},
  {Ukraine, "Ukraine", "UKR"},
  {United_States_Minor_Outlying_Islands, "United States Minor Outlying Islands", "UMI"},
  {Uruguay, "Uruguay", "URY"},
  {Uzbekistan, "Uzbekistan", "UZB"},
  {Vatican_City, "Vatican City", "VAT"},
  {Saint_Vincent_and_the_Grenadines, "Saint Vincent and the Grenadines", "VCT"},
  {Venezuela, "Venezuela", "VEN"},
  {British_Virgin_Islands, "British Virgin Islands", "VGB"},
  {US_Virgin_Islands, "Virgin Islands, U.S.", "VIR"},
  {Vietnam, "Vietnam", "VNM"},
  {Vanuatu, "Vanuatu", "VUT"},
  {Wallis_and_Futuna, "Wallis and Futuna", "WLF"},
  {Samoa, "Samoa", "WSM"},
  {Yemen, "Yemen", "YEM"},
  {South_Africa, "South Africa", "ZAF"},
  {Zambia, "Zambia", "ZMB"},
  {Zimbabwe, "Zimbabwe", "ZWE"},
  {Akrotiri_and_Dhekelia, "Akrotiri and Dhekelia", "---"},
  {Caspian_Sea, "Caspian Sea", "---"},
  {Clipperton_Island, "Clipperton Island", "CPT"},
  {Kosovo, "Kosovo", "XKX"},
  {Northern_Cyprus, "Northern Cyprus", "---"},
  {Paracel_Islands, "Paracel Islands", "---"},
  {Spratly_Islands, "Spratly Islands", "---"},
  {No_land, "No land", "---"},
  {Australia, "Australia", "AUS"},
  {Brazil, "Brazil", "BRA"},
  {Canada, "Canada", "CAN"},
  {China, "China", "CHN"},
  {India, "India", "IND"},
  {Russia, "Russia", "RUS"},
  {United_States, "United States", "USA"}
};

static int findcountryname(const char *name,
                           const Countryname countryname[],
                           int ncountries)
{
  int i;
  for(i=0;i<ncountries;i++)
    if(!strcmp(name,countryname[i].abbrev))
      return countryname[i].code;
  return NOT_FOUND;
} /* 'findcountryname' */

static int compare(const Countryname *a,const Countryname *b)
{
  return strcmp(a->name,b->name);
}

static Bool findcountry(const int country[],int n,int c)
{
  int i;
  for(i=0;i<n;i++)
    if(country[i]==c)
      return TRUE;
  return FALSE;
} /* of 'findcountry' */

int main(int argc,char **argv)
{
  FILE *file,*grid,*out;
  int i,*country,n,version,country_version;
  char *endptr;
  Intcoord coord;
  Header header,gridheader,outheader;
  String headername;
  int code;
  Bool swap_country,swap_grid,rc,isregion;
  float fcoord[2];
  double dcoord[2];
  outheader.nyear=1;
  outheader.nstep=1;
  outheader.timestep=1;
  outheader.firstcell=0;
  outheader.order=0;
  outheader.ncell=0;
  outheader.nbands=2;
  version=country_version=READ_VERSION;
  if(argc>1 && !strcmp(argv[1],"-longheader"))
  {
    version=country_version=2;
    argc--;
    argv++;
  }
  if(argc>1 && !strcmp(argv[1],"-list"))
  {
    puts("List of country codes:\nCode Name");
    qsort(countrynames,NCOUNTRY,sizeof(Countryname),
          (int (*)(const void *,const void *))compare);
    for(i=0;i<NCOUNTRY;i++)
      printf("%s  %s\n",countrynames[i].abbrev,countrynames[i].name);
    return EXIT_SUCCESS;
  }
  if(argc<5)
  {
    fprintf(stderr,"Argument(s) missing.\n"
            USAGE,argv[0]);
    return EXIT_FAILURE;
  }
  file=fopen(argv[1],"rb");
  if(file==NULL)
  {
    fprintf(stderr,"Error opening '%s': %s.\n",argv[1],strerror(errno));
    return EXIT_FAILURE;
  }
  if(freadanyheader(file,&header,&swap_country,headername,&country_version,TRUE))
  {
    fprintf(stderr,"Error reading header of '%s'.\n",argv[1]);
    return EXIT_FAILURE;
  }
  if(header.nbands==1)
    isregion=FALSE;
  else if(header.nbands!=2)
  {
    fprintf(stderr,"Invalid number of bands=%d in `%s', must be 1 or 2.\n",
            header.nbands,argv[1]);
    return EXIT_FAILURE;
  }
  else
   isregion=TRUE;
  grid=fopen(argv[2],"rb");
  if(grid==NULL)
  {
    fprintf(stderr,"Error opening '%s': %s.\n",argv[2],strerror(errno));
    return EXIT_FAILURE;
  }
  if(freadheader(grid,&gridheader,&swap_grid,LPJGRID_HEADER,&version,TRUE))
  {
    fprintf(stderr,"Error reading header of '%s'.\n",argv[2]);
    return EXIT_FAILURE;
  }
  outheader.firstyear=gridheader.firstyear;
  outheader.datatype=gridheader.datatype;
  if(version>=2)
  {
    outheader.cellsize_lon=gridheader.cellsize_lon;
    outheader.cellsize_lat=gridheader.cellsize_lat;
    outheader.scalar=gridheader.scalar;
  }
  n=argc-4;
  country=newvec(int,n);
  check(country);
  for(i=0;i<n;i++)
  {
    country[i]=strtol(argv[4+i],&endptr,10);
    if(*endptr!='\0')
    {
      /* argument is not a number */
      country[i]=findcountryname(argv[4+i],countrynames,NCOUNTRY);
      if(country[i]==NOT_FOUND)
      {
        fprintf(stderr,"Invalid number/name '%s' for country.\n",argv[4+i]);
        return EXIT_FAILURE;
      }
    }
  }
  out=fopen(argv[3],"wb");
  if(out==NULL)
  {
    fprintf(stderr,"Error creating '%s': %s.\n",argv[3],strerror(errno));
    return EXIT_FAILURE;
  }
  fwriteheader(out,&outheader,LPJGRID_HEADER,version);
  for(i=0;i<header.ncell;i++)
  {
    if(readcountrycode(file,&code,header.datatype,isregion,swap_country))
    {
      fprintf(stderr,"Error reading country code at %d.\n",i+1);
      return EXIT_FAILURE;
    }
    switch(gridheader.datatype)
    {
      case LPJ_SHORT:
        rc=readintcoord(grid,&coord,swap_grid);
        break;
      case LPJ_FLOAT:
        rc=freadfloat(fcoord,2,swap_grid,grid)!=2;
        break;
      case LPJ_DOUBLE:
        rc=freaddouble(dcoord,2,swap_grid,grid)!=2;
        break;
      default:
        fprintf(stderr,"Invalid datatype %d in '%s'.\n",gridheader.datatype,argv[2]);
        return EXIT_FAILURE;
    }
    if(rc)
    {
      fprintf(stderr,"Error reading coordinate at %d.\n",i+1);
      return EXIT_FAILURE;
    }
    if(findcountry(country,n,code))
    {
      switch(gridheader.datatype)
      {
        case LPJ_SHORT:
          rc=fwrite(&coord,sizeof(coord),1,out)!=1;
          break;
        case LPJ_FLOAT:
          rc=fwrite(fcoord,sizeof(float),2,out)!=2;
          break;
        case LPJ_DOUBLE:
          rc=fwrite(dcoord,sizeof(double),2,out)!=2;
          break;
        default:
          break;
      }
      if(rc)
      {
        fprintf(stderr,"Error writing coordinate at %d.\n",i+1);
        return EXIT_FAILURE;
      }
      outheader.ncell++;
    }
  }
  fclose(file);
  fclose(grid);
  rewind(out);
  fwriteheader(out,&outheader,LPJGRID_HEADER,version);
  fclose(out);
  if(header.ncell)
    printf("Number of cells: %d\n",outheader.ncell);
  else
    fputs("Warning: no cells found.\n",stderr);
  return EXIT_SUCCESS;
} /* of 'main' */
