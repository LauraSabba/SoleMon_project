
# Folder structure

If you want to use the procedure outside of the github folder
(recommended), you should replicate the same folder tree. In figure 1 is
provided an example of the structure that should be replicated. **Please
ensure that all the folders displayed in figure 1 exists**.

  

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/folder_str.png"
alt="Figure 1: folders tree structure" />
<figcaption aria-hidden="true">Figure 1: folders tree
structure</figcaption>
</figure>

# Access file creation and structure

## Starting a new access file

Data collected onboard are stored in access files, located in the
“OnBoard/assess” subfolder. Ideally, a maximum of 20 hauls should be
stored in a .accdb file.

To start a new file, copy the template version “Maschera inserimento
SOLEMON_template.accdb”, paste it in the same folder and rename it as
“Maschera inserimento SOLEMON_YEAR_N.accdb”, where YEAR is the reporting
year and N is a progressive number. Example: “**Maschera inserimento
SOLEMON_2024_1.accdb**” . *It does not matter if parts of the same haul
are in different files, this is handled in post-processing*.

## Starting a new haul

Each access file contains a template sheet, called “cala_template”. To
start a new haul, copy the template sheet, paste it and rename it as
“cala_x”, where *x* is the haul number. Example: **cala_1**; cala_7bis.
*It does not matter if parts of the same haul are in different files,
this is handled in post-processing*.

### Haul sheet Structure

The structure of the sheets, in 2024, is according to Figure1. In
detail:  

- **gear**: gear code. Accepts A and D.
- **species_name**: solemon code
- **length_mm**: length of the individual, in millimeters
- **weight_g**: weight of the individuals, in grams. This also host the
  cumulative weight for some species, refer to “deal with cumulative
  data” section for further details.
- **Sex**: required only for target species, accepts F, M, I. Leave
  empty if sex data not required
- **Mat**: required only for target species, just specify the stage (1,
  2, etc..). Exceptions for crustaceans: refer to section “cr”
- **id_specimen**: fill for specimens for which detailed samples
  (otolit, genetic etc.) were taken. It accepts numbers and letters.
  When filled a serial “fishID” number would be generated.
- **length_field2**: disc length for elasmobranchs; free spot for other
  species
- **length_field3**: disc width for elasmobranchs; free spot for other
  species
- **number_field1**: placeholder for cumulative number and subsamples.
  The compilation depends on the case
- **number_field1**: placeholder for cumulative number and subsamples.
  The compilation depends on the case
- **kg_field1**: weight in kilograms, placeholder for cumulative number
  and subsamples. The compilation depends on the case
- **kg_field2**: weight in kilograms, placeholder for cumulative number
  and subsamples. The compilation depends on the case
- **kg_field2**: weight in kilograms, placeholder for cumulative number
  and subsamples. The compilation depends on the case
- **type_subsample**: accepts values “haul”, “species”, “multi”,
  “other”. The compilation depends on the case.
- **Notes**: accept notes of any kind.

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/header_acces.png"
alt="Figure 2: Example of the structure of the access file" />
<figcaption aria-hidden="true">Figure 2: Example of the structure of the
access file</figcaption>
</figure>

### .accdb auto compilation general details

Auto compilation applies to some column of the access file. This mean
that in the processing of the table, to empty cells is assigned the
first available data that is found in the previous rows. When collecting
data, for these column you need to specify the value just for the first
observation, then you should fill the value again only when this change.
These columns are:

- **gear**
- **species_name**
- **sex** (only for crustaceans)

# Collecting data

## Standard procedure

Standard procedure applies when all the individuals for a given species
are collected and reported

### Target species

When reporting data for target species you should use only the first 7
columns. Of these, the columns *gear* and *species_name* autocompiles
according to the first record inserted. This mean that you should write
in these columns only the first time you report an observation (i.e.:
first record of a gear, first record of a species). The column
*id_specimen* serves to store any kind of individual id (e.g.: otholits
code, genetic samples).

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/target.png"
style="display: block; margin: 0 auto" width="600" height="100"
alt="Figure 3: example of compilation for target species species" />
<figcaption aria-hidden="true">Figure 3: example of compilation for
target species species</figcaption>
</figure>

### Elasmobranchs

When reporting data for elasmobranchs you should use only the first 9
columns. Of these, the columns *gear* and *species_name* autocompiles
according to the first record inserted. This mean that you should write
in these columns only the first time you report an observation (i.e.:
first record of a gear, first record of a species). The column
*id_specimen* serves to store any kind of individual id (e.g.: otholits
code, genetic samples). Refer to the Figure 4 for reporting the 3 length
measures

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/elas.png"
style="display: block; margin: 0 auto" width="500" height="400"
alt="Figure 4: example of compilation for elasmobranchs species" />
<figcaption aria-hidden="true">Figure 4: example of compilation for
elasmobranchs species</figcaption>
</figure>

### Other commercial species

When reporting data for other commercial species you should use only the
first 4 columns. Of these, the columns *gear* and *species_name*
autocompiles according to the first record inserted. This mean that you
should write in these columns only the first time you report an
observation (i.e.: first record of a gear, first record of a species).

For this species category are recorded individual length and cumulative
weight. The cumulative weight should be reported in the last record, as
in Figure 5.

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/other.png"
style="display: block; margin: 0 auto" width="500" height="100"
alt="Figure 5: example of compilation for other commercial species" />
<figcaption aria-hidden="true">Figure 5: example of compilation for
other commercial species</figcaption>
</figure>

### Shellfishes (when all the individuals are sampled)

When reporting data for shellfishes (MUREBRA, HEXATRU, OSTREDU) you
should use the columns ax reported in Figure 6. Of these, the columns
*gear* and *species_name* autocompiles according to the first record
inserted. This mean that you should write in these columns only the
first time you report an observation (i.e.: first record of a gear,
first record of a species).

NB: in case of subsamples, refer to the dedicated section

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/shell.png"
style="display: block; margin: 0 auto" width="1200" height="50"
alt="Figure 6: example of compilation for shellfishes species" />
<figcaption aria-hidden="true">Figure 6: example of compilation for
shellfishes species</figcaption>
</figure>

## Subsamples

Subsamples are taken when the amount of individuals in a given species
is too high to be processed. Considering the onboard practice used in
the solemon survey there are three cases of subsamples happening. Figure
7 reports the steps that are done according to the three cases. Each
case is treated differently depending on the type of species, and the
data treatment is explained in the dedicated sections.

In terms of onboard procedure, the cases general refers to:

- Case 1: ALL the individuals from a haul are collected (sorted sample),
  then the sorted sample is subsampled and processed
- Case 2: an unsorted subsample is taken from the haul, then it is
  sorted and ALL the individuals from the sorted subsample are processed
- Case 3: it applies for those shellfish coming in huge quantity
  (MUREBRA, HEXATRU). The subsampling happens in multiple steps (refer
  to dedicated section)

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/subsample_general.png"
style="display: block; margin: 0 auto" width="900" height="500"
alt="Figure 7: type of subsamples" />
<figcaption aria-hidden="true">Figure 7: type of subsamples</figcaption>
</figure>

### Target species

A typical case for subsample of target species is AEQUOPE, which usually
occurs in large quantities and LFD is needed. In this case only a few
individuals are measured to obtain the length structure, then it is
needed to estimate total number (and sometimes total weight).
Individuals that are processed for length and weight can be treated as
any other target species. Regarding the other information needed to
raise the values, you need to create a new record for the species (and
gear) to store subsample data. There are two expected cases:

#### **case 1 - Target**

ALL individuals collected from the haul (sorted sample), then a
subsample is taken (sorted subsample) for individual processing. One
record is dedicated to the subsample details. The total weight of the
sorted sample (in kilograms) is reported in the in the *kg_field1*
field. The weight of the sorted subsample is reported in the
*kg_field2*. *type_subsample* is “species”. The individuals subsampled
are processed according to the standard procedure for target species,
creating a new record for each specimen.

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/target_sub1.png"
style="display: block; margin: 0 auto" width="900" height="400"
alt="Figure 8: subsamples, case 1 for target species species" />
<figcaption aria-hidden="true">Figure 8: subsamples, case 1 for target
species species</figcaption>
</figure>

<br> <br>

#### **case 2 - Target**

A subsample is taken from the haul (unsorted subsample), then ALL the
individuals in the subsample are collected (sorted subsample) for
individual processing. One record is dedicated to the subsample details.
The weight of the haul is reported in *kg_field1*. The weigth of the
unsorted subsample is reported in *kg_field2*. The weight of the sorted
subsample is reported in *kg_field3*. *type_subsample* is “haul”. The
individuals subsampled are processed according to the standard procedure
for target species, creating a new record for each specimen.

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/target_sub2.png"
style="display: block; margin: 0 auto" width="900" height="400"
alt="Figure 9: subsamples, case 2 for target species species" />
<figcaption aria-hidden="true">Figure 9: subsamples, case 2 for target
species species</figcaption>
</figure>

<br> <br>

### Non-target species (MUREBRA, HEXATRU, etc..)

A typical case for subsample of non-target species is MUREBRA. When it
occurrs in large aggregations, total number (and sometimes total weight)
are estimated from subsamples. To store information needed to raise the
values, there are two expected cases:

<br>

#### **Case 1 - Non target**

ALL individuals collected from the haul (sorted sample), then a
subsample is taken (sorted subsample) for estimating total number. One
record is dedicated to the subsample details. The total weight of the
sorted sample (in kilograms) is reported in the in the *kg_field1*
field. The weight of the sorted subsample is reported in the
*kg_field2*. The number of individual in the sorted subsample is
reported in the *number_field1*. *type_subsample* is “species”.

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/nontarget_sub1.png"
style="display: block; margin: 0 auto" width="900" height="400"
alt="Figure 10: subsamples, case 1 for non target species species" />
<figcaption aria-hidden="true">Figure 10: subsamples, case 1 for non
target species species</figcaption>
</figure>

<br> <br>

#### **Case 2 - Non target**

A subsample is taken from the haul (unsorted subsample), then ALL the
individuals in the subsample are collected (sorted subsample) for
estimating total number. One record is dedicated to the subsample
details. The weight of the haul is reported in *kg_field1*. The weigth
of the unsorted subsample is reported in *kg_field2*. The weight of the
sorted subsample is reported in *kg_field3*. The number of individual in
the sorted subsample is reported in the *number_field1*.
*type_subsample* is “haul”.

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/nontarget_sub2.png"
style="display: block; margin: 0 auto" width="900" height="400"
alt="Figure 11: subsamples, case 2 for non target species species" />
<figcaption aria-hidden="true">Figure 11: subsamples, case 2 for non
target species species</figcaption>
</figure>

<br> <br>

#### **Case 3 - Non target**

it has been only used for MUREBRA and HEXATRU in cases of exceptional
catches. ALL individuals of these two species are collected from the
haul (partially sorted sample), then a subsample is taken from the
partially sorted sample (partially sorted subsample). The partially
sorted subsample is sorted, aka it is divede into species (sorted
subsample). The sorted subsample of each species is subsampled againg
(sorted sub subsample) for estimating total number. One record is
dedicated to the subsample details. The weight of the partially sorted
sample is reported in *kg_field1*. The weigth of the sorted subsample
for each species is reported in *kg_field2*. The weight of the sorted
sub-subsample is reported in *kg_field3*. The number of individual in
the sorted sub-subsample is reported in the *number_field1*.
*type_subsample* is “multi”.

<figure>
<img
src="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/Solemon/github/SoleMon_project/OnBoard/data/other/nontarget_sub3.png"
style="display: block; margin: 0 auto" width="900" height="400"
alt="Figure 12: subsamples, case 3 for non target species species" />
<figcaption aria-hidden="true">Figure 12: subsamples, case 3 for non
target species species</figcaption>
</figure>

<br> <br>

# Processing data

Data are processed in R. The code runs on R 4.0.5 32 bit version due to
compatibility issue with RODBC package needed to read from .accdb. It is
good practice to process the data relatively frequently in order to
create backup files and to produce some checks that can be used to
adjust data.

## process haul data

Data stored in the .accdb file are retrieved and handled by R scripts,
located in the “R” folder. The required script is `workflow_access_v0`.

### process single hauls

To process single hauls, run the lines up to 15 and then please set the
following parameters:

- **haul**: number of the haul corresponding to the name given in the
  access file, but without the extension ’‘. Example: ’cala22’ in access
  is ‘22’ here, ‘cala45bis’ is ‘45bis’
- **db**: suffix given to the acces file. E.g.: if access file name is
  *‘Maschera inserimento SOLEMON_test.accdb’*, type just *‘test’* here.
- **updateID**: this control if update (‘Y’) or not (‘N’) the serial
  number ‘fishID’. Recommended to type ‘Y’ only when batch processing
  all the hauls.
- **area_sepia**: old command
- **year**: year of the survey
- **area**: stratum of the haul

<!-- -->

    # set parameters
    haul=22 
    db='test'
    updateID='N'
    area_sepia='D'
    year=2022
    area='ITA17'

  

When parameters are set, the data processing is pre-defined: you just
have to run it without changing the parameters. The first step is the
`function1`, which scope is to format the access table according to
output standards. There is no need to see the output of this function.
Just run it.

    # function1 extract data from access db and format them
    hauldata=function1(haul=haul, 
                       db=db,
                       year=year)# extract and format data

  

When the function1 is done, you can proceed with the `function2`, which
scope is to perform some checks. Checks done are plots that would be
saved under the path ‘output/checks’

    # function 2: perform checks
    function2(xdat=hauldata, 
              haul=haul)

  

`function3` scope is to format data according to trust format. Excel
sheets are saved under the path ‘output/trust’.

    # function 3: format data to trust format
    trustdat=function3(xdat=hauldata[[1]], 
                      haul=haul, 
                      year = year, 
                      weight_not_target = hauldata[[2]],  
                      subsamples_target=hauldata[[3]],
                      catch_sample_disattivati = catch_sample_disattivati) # function 2
                      
                      

  

`function4` creates a pdf report and save it under the path
‘output/pdf’.

    # function4: save PDF
    function4(trustdat = trustdat, 
              year=year,
              area = area,
              haul=haul)

### loop to process more hauls

To process more than one haul, you should care to properly fill in the
‘haul_order’ excel sheet (see input data section). After having loaded
the haul summary, just run the loop represented below.

    haul_summary=read_excel("data/haul_order.xlsx")
    haul_summary=haul_summary[1:5,]

    for(xhaul in 1:nrow(haul_summary)){
      
      
      # loop parameters
      haul=haul_summary[xhaul,]$haul
      db=haul_summary[xhaul,]$DB
      area=haul_summary[xhaul,]$country
      
      cat('processing haul no.', haul, '(', xhaul,'/', nrow(haul_summary),')' )
      
      # function1 extract data from access db and format them
      hauldata=function1(haul=haul, 
                         db=db,
                         year=year)# extract and format data
      
      # function 2: perform checks
      function2(xdat=hauldata, 
                haul=haul)
      
      # function 3: format data to trust format
      trustdat=function3(xdat=hauldata[[1]], 
                         haul=haul, 
                         year = year, 
                         weight_not_target = hauldata[[2]],  
                         subsamples_target=hauldata[[3]],
                         catch_sample_disattivati = catch_sample_disattivati) # function 2
      
      # function4: save PDF
      function4(trustdat = trustdat, 
                year=year,
                area=area,
                haul=haul)


    }

# Data files associated to the folder

### target_species

Contains the species that are target for the survey (individual length
and weight) and the molluscs for which only total weight and total
number are needed.

- **Species**: species solemon code
- **target**: indicates which type of target is the species. 1 is for
  species that requires individual length and weight; 2 is for molluscs
  that requires total weight and total number

  
  

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>
List of target species. Target = 1 are proper target; Target = 2 are
shellfishes on which to collect number and weight
</caption>
<thead>
<tr>
<th style="text-align:left;">
Species
</th>
<th style="text-align:left;">
Species_code
</th>
<th style="text-align:right;">
target
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Aequipecten opercularis
</td>
<td style="text-align:left;">
AEQUOPE
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Callinectes sapidus
</td>
<td style="text-align:left;">
CALLSAP
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Flexopecten glaber proteus (= Chlamys proteus)
</td>
<td style="text-align:left;">
CHLAPRO
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Chlamys varia
</td>
<td style="text-align:left;">
CHLAVAR
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Melicertus kerathurus
</td>
<td style="text-align:left;">
MELIKER
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Merluccius merluccius
</td>
<td style="text-align:left;">
MERLMER
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Mullus barbatus
</td>
<td style="text-align:left;">
MULLBAR
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Nephrops norvegicus
</td>
<td style="text-align:left;">
NEPRNOR
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Parapenaeus longirostris
</td>
<td style="text-align:left;">
PAPELON
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Pecten jacobaeus
</td>
<td style="text-align:left;">
PECTJAC
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Platichthys flesus
</td>
<td style="text-align:left;">
PLATFLE
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Psetta maxima
</td>
<td style="text-align:left;">
PSETMAX
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Raja asterias
</td>
<td style="text-align:left;">
RAJAAST
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Raja clavata
</td>
<td style="text-align:left;">
RAJACLA
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Raja miraletus
</td>
<td style="text-align:left;">
RAJAMIR
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Scophthalmus rhombus
</td>
<td style="text-align:left;">
SCOHRHO
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Scyliorhinus canicula
</td>
<td style="text-align:left;">
SCYOCAN
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Scyliorhinus stellaris
</td>
<td style="text-align:left;">
SCYOSTE
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sepia officinalis
</td>
<td style="text-align:left;">
SEPIOFF
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Solea aegyptiaca
</td>
<td style="text-align:left;">
SOLEAEG
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Solea solea
</td>
<td style="text-align:left;">
SOLEVUL
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Squilla mantis
</td>
<td style="text-align:left;">
SQUIMAN
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Torpedo marmorata
</td>
<td style="text-align:left;">
TORPMAR
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Crassostrea gigas
</td>
<td style="text-align:left;">
CRASGIG
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Galeodea echinophora
</td>
<td style="text-align:left;">
GALEECH
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Hexaplex trunculus
</td>
<td style="text-align:left;">
HEXATRU
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Bolinus brandaris
</td>
<td style="text-align:left;">
MUREBRA
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Mytilus galloprovincialis
</td>
<td style="text-align:left;">
MYTGALL
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Natica stercusmuscarum
</td>
<td style="text-align:left;">
NATISTE
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Ostrea edulis
</td>
<td style="text-align:left;">
OSTREDU
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Rapana venosa
</td>
<td style="text-align:left;">
RAPAVEN
</td>
<td style="text-align:right;">
2
</td>
</tr>
</tbody>
</table>

### catch_sample_disattivati

This file controls formatting of trust templates. It indicates which
samples (Station, Gear and SpecCode) should be indicated as “InUse” =
FALSE in the catch sample files used as input data in trust.

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
Station
</th>
<th style="text-align:left;">
Gear
</th>
<th style="text-align:left;">
SpecCode
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
D
</td>
<td style="text-align:left;">
AEQUOPE
</td>
</tr>
</tbody>
</table>

  
  

### fishID

Store the *updated* serial number used to identify specimens for which
detailed samples (otolit, genetic etc.) were taken. This number should
refere to the *last* ID assigned to a specimen. The use of this file is
controlled by the `updateID` parameter in the workflow_access_v0 file:
if `updateID` is set equal to Y, the fishID file is used to assign IDs
(when requested) and it is then updated. The columns refers to:

- **type**: not useful (to be removed)
- **code**: this is the alphanumeric part of the code to be assigned to
  specimens
- **fishID**: serial number referring to the *last* individual contained
  in the past records
- **haul**: not useful (to be removed)
- **species**: indicates which species belong to the code category. If
  need to add species, use : as separator and do not insert spaces.

  
  

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>
Example of the fishID file structure
</caption>
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
code
</th>
<th style="text-align:right;">
fishID
</th>
<th style="text-align:left;">
haul
</th>
<th style="text-align:left;">
species
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ELAS
</td>
<td style="text-align:left;">
Elas
</td>
<td style="text-align:right;">
202
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
RAJAAST:RAJACLA:RAJAMIR:TORPMAR
</td>
</tr>
<tr>
<td style="text-align:left;">
solea
</td>
<td style="text-align:left;">
SS
</td>
<td style="text-align:right;">
7061
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
SOLEVUL:SOLEAEG
</td>
</tr>
<tr>
<td style="text-align:left;">
RHO
</td>
<td style="text-align:left;">
SR
</td>
<td style="text-align:right;">
447
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
SCOHRHO
</td>
</tr>
</tbody>
</table>

  
  

### haul_order

Store the information associated with hauls. This file is used (1) when
the data workflow is applied in loop; (2) by the minilog script. The
columns refers to:

- **day**: day of the haul (yyyy-mm-dd)
- **haul**: number of the haul as from solemon protocol
- **id**: progressive number of the haul
- **note**: this space serve to write any kind of note, it is ignored by
  the code
- **inizio**: time of setting the net (hh:mm:ss)
- **fine**: time of hauling the net (hh:mm:ss)
- **verifica_shell**: this space serve to write additional notes, it is
  ignored by the code
- **DB**: indicates the name of the access database where the haul was
  recorder. Do not include ‘.accdb’ extension.
- **country**: indicates in which stratum the haul was performed.
  Available strata are ‘HRV’, ‘ITA17’, ‘SVN’

  
  

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>
Example of the haul_order file structure
</caption>
<thead>
<tr>
<th style="text-align:left;">
day
</th>
<th style="text-align:right;">
haul
</th>
<th style="text-align:right;">
id
</th>
<th style="text-align:left;">
note
</th>
<th style="text-align:left;">
inizio
</th>
<th style="text-align:left;">
fine
</th>
<th style="text-align:right;">
valid
</th>
<th style="text-align:left;">
DB
</th>
<th style="text-align:left;">
country
</th>
<th style="text-align:left;">
peso_rapido_A
</th>
<th style="text-align:left;">
peso_rapido_D
</th>
<th style="text-align:left;">
peso_subcampione_a
</th>
<th style="text-align:left;">
peso_subcampione_D
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
2022-11-24
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
2022_1
</td>
<td style="text-align:left;">
ITA17
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
2022-11-24
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
2022_1
</td>
<td style="text-align:left;">
ITA17
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
2022-11-24
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
2022_1
</td>
<td style="text-align:left;">
ITA17
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
</tbody>
</table>

  
  

### lw_pars

Store the length-weigth parameters of target species. This file is used
to reconstruct length (or weight) when it is not available in the
recorded data. Example: shrimps where missi a part of the tail but have
the head intact are suitable only for length measurement; fishes that
were spoiled by the gear may be ok for weight but not measurables for
length. The columns refers to:

- **species_name**: indicates which species belong to the code category
  (solemon code)
- **sex**: indicates the sex in case when sexual dimorphism in the
  growth is relevant. If dimorphism is not relevant, put combined
  parameter and write NA in this field
- **a**: parameter “a” for the length-weight relationship
- **b**: parameter “b” for the length-weight relationship
- **source**: source from where a and b were retrieved

  
  

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>
Example of the lw_pars file structure
</caption>
<thead>
<tr>
<th style="text-align:left;">
species_name
</th>
<th style="text-align:left;">
sex
</th>
<th style="text-align:right;">
a
</th>
<th style="text-align:right;">
b
</th>
<th style="text-align:left;">
source
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
SOLEVUL
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.00460
</td>
<td style="text-align:right;">
3.11
</td>
<td style="text-align:left;">
benchmark_assessment
</td>
</tr>
<tr>
<td style="text-align:left;">
MULLBAR
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.00871
</td>
<td style="text-align:right;">
3.09
</td>
<td style="text-align:left;">
fishbase
</td>
</tr>
<tr>
<td style="text-align:left;">
RAJACLA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.00269
</td>
<td style="text-align:right;">
3.23
</td>
<td style="text-align:left;">
fishbase
</td>
</tr>
</tbody>
</table>

  
  

### maturity_stages

Store the code of the maturity scales for target species. This file is
used to format input files for trust.

- **SPECIES**: indicates which species belong to the code category
  (solemon code)
- **SEX**: indicates the sex
- **SCALE**: include the alphanumeric part of the maturity scale code

  
  

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>
Example of the maturity_stages file structure
</caption>
<thead>
<tr>
<th style="text-align:left;">
SPECIES
</th>
<th style="text-align:left;">
SEX
</th>
<th style="text-align:left;">
SCALE
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
MELIKER
</td>
<td style="text-align:left;">
F
</td>
<td style="text-align:left;">
MEDPF
</td>
</tr>
<tr>
<td style="text-align:left;">
MERLMER
</td>
<td style="text-align:left;">
F
</td>
<td style="text-align:left;">
MEDFI
</td>
</tr>
<tr>
<td style="text-align:left;">
MERLMER
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:left;">
MEDFI
</td>
</tr>
</tbody>
</table>

  
  

### solemon_TB

This is just the TB file updated to 2021, as stored in trust. It serves
to perform some checks and it should not be modified. Preview not shown.

  
  

### species_list

This file is downloaded from the trust database and not modified. It
contains the species list. Last update xxx. If need to modify, please
*do it in thrust and then download the excel again!*.

- **Species**: species scientific name
- **Medits**: species solemon code
- **Sp_Subcat**: commercial category of the species
- **Lan_Class**: not relevant
- **Sp_Subcat**: not relevant

  
  

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>
Example of the species_list file structure
</caption>
<thead>
<tr>
<th style="text-align:left;">
Species
</th>
<th style="text-align:left;">
Medits
</th>
<th style="text-align:left;">
Sp_Subcat
</th>
<th style="text-align:left;">
Len_Class
</th>
<th style="text-align:left;">
Notes
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Aaptos aaptos
</td>
<td style="text-align:left;">
AAPTAAP
</td>
<td style="text-align:left;">
E
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Abra prismatica
</td>
<td style="text-align:left;">
ABRAPRI
</td>
<td style="text-align:left;">
E
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Abra sp
</td>
<td style="text-align:left;">
ABRASPP
</td>
<td style="text-align:left;">
E
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
</tbody>
</table>

  
  
