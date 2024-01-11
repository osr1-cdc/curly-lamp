# Tableau Util
This is a little set of code to access the CDC EDAV Tableau Server programatically. It is a long time coming. Tableau now has documentation on how to actually use the API here: https://help.tableau.com/current/api/rest_api/en-us/REST/rest_api.htm

The main thing is that you need to have a file called `signing.xml` in the same folder as the python script here. The contents of this file will include YOUR CDC username and password, so please be careful NOT to post it back here. I have added tableau_util/signin.xml to .gitignore just in case. Here is what should be in the signin.xml file:

```xml
<tsRequest>
 <credentials name="cdcusername" password="cdcpassword" >
   <site contentUrl="" />
 </credentials>
</tsRequest>
```

Importantly, each site, folder, workbook, and view can and must be referred to by a UUID (Tableau also calls these LUIDs). In this script, there are functions to find the correct workbook based on name. You should run by finding the exact name of the workbook on https://tableau.edav.cdc.gov/#/projects/36 then use it when calling this program:

```bash
python pull_tableau_data.py Variant_Proportions_Plus_Nowcasting_240103
```


It will automatically try to infer the dates you want for latest Nowcast and Weighted data, but you cal aso supply those as additional positional arguments:


```bash
python pull_tableau_data.py Variant_Proportions_Plus_Nowcasting_240103 2023-12-24 2024-01-20
```