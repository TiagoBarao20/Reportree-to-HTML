#!/usr/bin/env  python3

"""
This script was made with the intent to improve the user's visualization and accessbility when interpreting the reportree outputs 

By Tiago Barão
Github:https://github.com/TiagoBarao20
"""

import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import os

############# Functions #############
"this function gets the user's files of interest as input and stores them in a list"
def get_file_paths():
    file_paths = []
    while True:
        file_path = input("Please enter your files(or 'done' to finish): ")

        if file_path.lower() == "done":
            break
    
        if os.path.isfile(file_path):
            file_paths.append(file_path)
    return file_paths

def info_log(infile):
    """ this funcion retrieves info(date, run name, comand used and version) from the Reportree log file ("Lm.log")"""

    with open(infile, "r") as fp:
        line_numbers = [3, 5, 6]
        lines =  []
        for i, line in enumerate(fp):
            if i in line_numbers:
                lines.append(line.strip())
            elif i > 6:
                break
    version = lines[0]
    Command = lines[1]
    Date = lines[2]
    versionbreak = version.split()
    finalversion = versionbreak[1]
    datebreak = Date.split()
    finaldate = datebreak[1]
    run_name = r"-out\s+(.*?)\s+--analysis"
    output_path = re.search(run_name, Command).group(1)

    return Command, finalversion, finaldate, output_path

def get_header(Command, finalversion, finaldate, output_path):
    "This function gets as input the information we retrieved from the Reportree log file and uses it to make the head of the HTML file using"
    "the run name and also uses a png file with the Reportree logo"

    header = '''<!DOCTYPE html>
    <html lang="en">
    <head>
    <title>Reportree run: ''' + output_path + '''</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
    /* Style the body */
    body {
    font-family: Arial;
    margin: 0;
    }

    /* Header/Logo Title */
    .header {
    padding: 60px;
    text-align: center;
    background: #ffffff;
    color: rgb(0, 0, 0);
    font-size: 30px;
    }

    .div-img {
    height: 250px;
    width: 400px;
    background-size: cover;
    background-image: url("https://github.com/insapathogenomics/ReporTree/blob/main/reportree_logo.png?raw=true");
    background-position: center center;
    background-repeat: no-repeat;
    /* Center the div horizontally and vertically */
    margin: 0 auto;
    }

    /* Page Content */
    .content {
    padding: 20px;
    }
    </style>
    </head>
    <body>

    <div class="header">
    <h1></h1>
    <div class="div-img"></div>
    <p>Run name ''' + output_path + '''</p>
    </div>
    <div class="content">
    <h1>Parametros</h1>
    <p>Version: ''' + "".join(finalversion) + '''</p>
    <p>Command used: ''' + "".join(Command) + '''</p>
    <p>Run Date: ''' + "".join(finaldate) + '''</p>
    </div>

    </body>
    </html>
    '''
    return header

def particoes_clusters(infile):
    """ This function returns the number of clusters and the Threshold from the Reportree partitions table"""
    tabela = pd.read_table(infile)
    threshold_list = []
    nclusters_list = []
    count_column=0

    for column_name in tabela.columns[1:]:
        unique_values = tabela[column_name].nunique()
        threshold_list.append(str(count_column))
        nclusters_list.append(str(unique_values))
        count_column += 1
    
    threshold = "[" + ",".join(threshold_list) + "]"
    nclusters = "[" + ",".join(nclusters_list) + "]"
    highest_y = nclusters_list[0]
    size_x = len(threshold_list)

    return threshold, nclusters , highest_y, size_x

def get_lineplot(threshold, nclusters , highest_y, size_x):
    "This function uses as imput the number of clusters and their respective Threshold value to use in the creation of a Line Graph plot in the HTML"
    plot = '''
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                  
    <div id="myPlot" style="width:100%;max-width:700px"></div>

    <script>
        
    const xArray = ''' + threshold + ''';
    const yArray =  ''' + nclusters + ''';


    const data = [{
        x: xArray,
        y: yArray,
        mode:"lines"
    }];


    const layout = {
        xaxis: {range: [0,'''+(str(size_x))+'''], title: "Threshold"},
        yaxis: {range: [0, '''+(str(highest_y))+'''], title: "Number of clusters"},  
        title: "Number of cluster per Threshold value"
    };


    Plotly.newPlot("myPlot", data, layout);
    </script>
    <body>'''
    return plot

def nomenclature_changes(infile):
    "This function retrieves info from the nomenclature changes file, filters the information based on when the change is either new or increase"
    "and stores it into a DataFrame to later be used in the HTML"
    tabela = pd.read_csv(infile, sep='\t', header=None)
    filtered_tabela = tabela[(tabela[3].str.lower().str.contains('increase', na=False)|tabela[1].str.lower().str.contains('new', na=False))]           
    new_column_names = ["Partition", " Old Cluster", "Old Cluster Lenght", "Nomenclature Changes", "Cluster "+ finaldate + "","Cluster "+ finaldate + " Lenght", "Increased by", "Samples Increased"]
    filtered_tabela.columns = new_column_names
    return filtered_tabela

def get_table(filtered_tabela):
    "This function creates the HTML table using the filtered nomenclture changes DataFrame"
    tabelahtml = filtered_tabela.to_html(index=False)
    return tabelahtml

def metadata(infile, filtered_tabela):
    "This function retrieves info from both the metadata file and the filtered by samples of interest Nomenclatures Changes file"
    tabela = filtered_tabela
    tabelaparticoes = pd.read_csv(infile, sep='\t')
    common_values = set(tabelaparticoes.columns).intersection(tabela['Partition'])
    todas_x = []
    todas_y = []
    for matching_value in common_values:
        filtered_df = tabelaparticoes.loc[tabelaparticoes[matching_value].fillna('').str.contains('cluster', case=False)]
        for col, row in tabela.iterrows():
            if row['Partition'] == matching_value:
                cluster_name = filtered_df[matching_value].values[0]
                cluster_length = int(row[tabela.columns[5]])
                if cluster_length > 1:
                    flt_tabelaparticoes = tabelaparticoes[tabelaparticoes[matching_value] == cluster_name]
                    Coords = []
                    Coordsy = []
                    Coords = np.array(flt_tabelaparticoes['country'].unique())
                    Coords_list = list(Coords)
                    Coordsy = flt_tabelaparticoes['country'].value_counts()
                    Coordsy_list = Coordsy.tolist()
                    str_x= str(Coords_list)
                    str_y=str(Coordsy_list)
                    todas_x.append(str_x)
                    todas_y.append(str_y)
    
    return todas_x, todas_y 

def piechart(todas_x, todas_y):
    "This function uses the list of coordenates obtained from the previous function to create the dynamic variables used to create the pie chart plot"

    names = []
    scriptplot = ""
    Classesrepeat = ""
    for i in range(0, len(todas_x)):
        name = '''#myPlot'''+str(i)
        names.append(name)
        namepie = ",".join(names)
        Classesrepeat += '''<div id="myPlot'''+str(i)+'''"></div>\n'''
        todas_x_val=todas_x[i]
        todas_y_val =todas_y[i]
        todas_x_str = str(todas_x_val)
        todas_y_str = str(todas_y_val)
        scriptplot += '''
        <script>

        const xArray'''+str(i)+''' = ''' + todas_x_str + ''';
        const yArray'''+str(i)+''' = ''' + todas_y_str + ''';

        const layout'''+str(i)+''' = {title:"Country destribuition"};

        const data'''+str(i)+''' = [{labels:xArray'''+str(i)+''', values:yArray'''+str(i)+''', type:"pie"}];

        Plotly.newPlot("myPlot'''+str(i)+'''", data'''+str(i)+''', layout'''+str(i)+''');
        </script> \n
        '''
    return scriptplot, namepie, Classesrepeat
    
def get_pie(Classesrepeat):
    "this function uses one of the dynamic variables obtained in the function prior to this one and creates static parts of the pie chart plot"
    sizepie = ''' { ''''''
      width: 100%;
      max-width: 700px;
      height: 300px;
      }
    </style>'''
    headpie = '''<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        .charts-container {
        display: flex;
        flex-wrap: wrap;
        }
    '''
    Classespie = '''<body>
    <div class="charts-container">\n
    '''+Classesrepeat+'''
    </div>
    '''
    return sizepie, headpie, Classespie


def CreateHTML(lista):
    """This function recieves as input the list of variables necessary to create the HTML file and creates it"""
    with open("Report.html", "w") as tst:
        tst.write("".join(lista))

def sample_table(infile):
    "This function recieves as input the samples of interest file, changes its column names and creates the varaible to generate the HTML table"
    tabela = pd.read_csv(infile, delimiter="\t")
    new_column_names = ["Samples of Interest", "Partition", "Cluster", "Nomenclature Changes", "Increased by","Cluster Lenght", "Samples","Country"," Number of Countries","Source", "Samples Increased"]
    tabela.columns = new_column_names
    tabelasample = tabela.to_html(index=False)
    return tabelasample

def get_tabsample(tabelasample):
    "This function recieves the variable from the samples of interest and creates the HTML variable to create the table"
    styled_table = f'''
    <style>
        table {{
            table-layout: fixed;
            width: 85%;
            border-collapse: collapse;
        }}

        th, td {{
            padding: 8px;
            text-align: left;
            border: 1px solid #ddd;
            word-wrap: break-word; /* Enable text wrapping */
        }}

        /* Adjust column width for 'samples' column */
        th:nth-child(7), td:nth-child(7) {{ /* '1' corresponds to the first column */
            width: 300px; /* Adjust the width as needed */
        }}
        th:nth-child(4), td:nth-child(4) {{ /* '1' corresponds to the first column */
            width: 150px; /* Adjust the width as needed */
        }}

        /* Add more nth-child rules for other columns if necessary */
        .table-container {{
        margin-top: 30px;
        }}
    </style>
    <div class ="table-container">
        <table>
        {tabelasample}
        </table>
    </div>
    '''
    return styled_table

def spacetable(styled_table):
    "This function is used to separate all the tables equally in the HTML file"
    html_tables = styled_table
    return html_tables

############# Pipeline #############
file_paths = get_file_paths()
print(file_paths)
variables_for_html = []

#Objective1: Obtain a header variable containing the reportree logo and the run name
#Input:Png file for the logo and Log file
#Output:header variable

if "Lm.log" in file_paths:
    log_file_index = file_paths.index("Lm.log")
    log_file_path = file_paths[log_file_index]
    Command, finalversion, finaldate, output_path = info_log(log_file_path)
    header = get_header(Command, finalversion, finaldate, output_path)
    variables_for_html.append(header)
else:
    print("Lm.log file not provided. Skipping the processing of 'Lm.log'.")

#Objective2: Obtain the Line Plot variable
#Input:Partitions File
#Output:Line Plot variable

if "Lm_partitions.tsv" in file_paths:
    tsv_file_index = file_paths.index("Lm_partitions.tsv")
    file_path2 = file_paths[tsv_file_index]
    threshold, nclusters , highest_y, size_x = particoes_clusters(file_path2)
    lineplot = get_lineplot(threshold, nclusters , highest_y, size_x)
    variables_for_html.append(lineplot)
else:
    print("Lm_partitions.tsv file not provided. Skipping the processing.")

#Objective3: Obtain a variable containing cluster nomenclature changes for the new samples(only for increase and new)
#Input:Nomenclature Changes File
#Output: Variable with cluster nomenclature changes
if "Lm_nomenclature_changes.tsv" in file_paths:
    tsv_file_index = file_paths.index("Lm_nomenclature_changes.tsv")
    file_path3 = file_paths[tsv_file_index]
    filtered_tabela = nomenclature_changes(file_path3)
    tabelahtml = get_table(filtered_tabela)
    variables_for_html.append(tabelahtml)
else:
    print("Lm_nomenclature_changes.tsv file not provided. Skipping the processing.")

#Objective4: Obtain the variable to recreate the Samples of Interest file as a table in the HTML File
#Input: Samples of Interest with partitions file
#Output:Variable to recreate the table

if "Lm_SAMPLES_OF_INTEREST_partitions_summary.tsv" in file_paths:
    tsv_file_index = file_paths.index("Lm_SAMPLES_OF_INTEREST_partitions_summary.tsv")
    file_path4 = file_paths[tsv_file_index]
    tabelasample = sample_table(file_path4)
    styled_table = get_tabsample(tabelasample)
    html_tables = spacetable(styled_table)
    variables_for_html.append(html_tables)
else:
    print("Lm_SAMPLES_OF_INTEREST_partitions_summary.tsv file not provided. Skipping the processing.")


#Objective5: Obtain the variables necessary for the pie chart plot
#Input:Metadata with partitions file and Filtered table from the last objective
#Output: Variable with Pie Chart instructions

if "Lm_metadata_w_partitions.tsv" in file_paths:
    meta_file_index = file_paths.index("Lm_metadata_w_partitions.tsv")
    infile = file_paths[meta_file_index]
    todas_x, todas_y = metadata(infile, filtered_tabela)
    scriptplot, namepie, Classesrepeat = piechart(todas_x, todas_y)
    sizepie, headpie, Classespie = get_pie(Classesrepeat)
    pie_chart_complete = headpie + namepie + sizepie + Classespie + scriptplot
    print (pie_chart_complete)
    variables_for_html.append(pie_chart_complete)
else:
    print("Lm_metadata_w_partitions.tsv file not provided. Skipping the processing.")

#Objective6: Create the HTML File
#Input: Every output variable from previous objectives
#Output: HTML File

CreateHTML(variables_for_html)