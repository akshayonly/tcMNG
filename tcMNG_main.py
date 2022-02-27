"""
Title: Topics Compare MeSH Network Graph (tcMNG)
Author: Akshay Shirsath
Logo Icon Source: flaticon.com
"""
############################################
############### Libraries ##################

import streamlit as st
import streamlit.components.v1 as components
from PIL import Image
from string import Template
import datetime
import os 

import CustomFunctions
from pyvis.network import Network

############################################
############### Page Body ##################

st.set_page_config(layout="wide")

path = "tcMNG_logo.png"

image = Image.open(path)
st.image(image, use_column_width=True)

# st.title('MeSH Knowledge Graph')

st.markdown(
    """
    Visualise & Compare Any Two Research (Scientific) Topics Based On MeSH Terms.
    """)

expander_bar = st.expander("What is tcMNG?")
expander_bar.markdown("""
Topics Comparison MeSH Network-Graph

This web application takes two query terms. With some parameters fetches their respective PubMed records, extracts MeSH terms, and then builds Network-graph. 

The idea is to explore diverse or distinct topics and see what similarities they can bring forth. 
The tool is not meant to replace the reading process, but to show the insightful connections between different articles, and then the researcher can invest his further time for reading only relevant ones.
""") 

expander_bar = st.expander("More Info")
expander_bar.markdown("""
* **version 1.0**
* **Python Libraries:** streamlit, Pyvis & Biopython.
* **Data source:** NCBI Entrez
* **Author:** Akshay Shirsath  
* **Check out Network-Graph for your own article** [Click-Here](https://samng.herokuapp.com/)
""")    

default_one = 'Cancer'
default_two = 'SARS-CoV-2'

st.subheader('Compare')
query_one = st.text_input('Topic I', default_one)

st.subheader('And')
query_two = st.text_input('Topic II', default_two)

############### Parameters ##################

st.subheader('Parameters')

date_min = st.radio('Publication Date', ['1 year', '5 years', '10 years'])

present_date = datetime.datetime.now()

default_max = (str(present_date.year)+"/"+str(present_date.month)+"/"+str(present_date.day))

if date_min == '1 year':
    default_min = (str(present_date.year-(1))+"/"+str(present_date.month)+"/"+str(present_date.day))
if date_min == '5 years':
    default_min = (str(present_date.year-(5))+"/"+str(present_date.month)+"/"+str(present_date.day))
if date_min == '10 years':
    default_min = (str(present_date.year-(10))+"/"+str(present_date.month)+"/"+str(present_date.day))

article_retrive = st.radio('Articles To Retrive (for each topics)', [3, 5, 10])

if st.button('Show Graph'):
    message_one = f"Searching PMIDs associated with Topic I ({query_one})..."
    with st.spinner(text=message_one):
        results_one = CustomFunctions.fetch_pmids(query_one, article_retrive)

    message_two = f"Searching PMIDs associated with Topic II ({query_two})..."
    with st.spinner(text=message_two):
        results_two = CustomFunctions.fetch_pmids(query_two, article_retrive)        

    with st.spinner(text=f"Processing Topic I ({query_one}) data..."):
        pmids_one = ",".join(results_one['IdList'])
        articles_data_one = CustomFunctions.fetch_articles(pmids_one, min_date=default_min, max_date=default_max)

    with st.spinner(text=f"Processing Topic II ({query_two}) data..."):
        pmids_two = ",".join(results_two['IdList'])
        articles_data_two = CustomFunctions.fetch_articles(pmids_two, min_date=default_min, max_date=default_max)   

    MKGraph = Network(height='700px', width='81%', bgcolor='#222222', font_color='#ecf0f1')
    
    with st.spinner(text=f"Searching Common Terms..."):
        common_terms = CustomFunctions.common_mesh_terms(articles_data_one, articles_data_two)
        
    with st.spinner(text=f"Generating Graph for Topic I ({query_one})..."):
        MKGraph = CustomFunctions.create_graph(articles_data_one, MKGraph, common_terms, color_type='FIRST')        

    with st.spinner(text=f"Generating Graph for Topic II ({query_two})..."):
        MKGraph = CustomFunctions.create_graph(articles_data_two, MKGraph, common_terms, color_type='SECOND')          

        MKGraph.set_options("""
        var options = {
        "edges": {
            "arrows": {
            "to": {
                "enabled": true,
                "scaleFactor": 0.5
            }
            },
            "color": {
            "inherit": true
            },
            "smooth": {
            "forceDirection": "none"
            }
        },
        "physics": {
            "barnesHut": {
            "gravitationalConstant": -17350,
            "springLength": 210,
            "springConstant": 0.055,
            "avoidOverlap": 0.53
            },
            "minVelocity": 0.75
        }
        }
        """)

        t = Template(
                """
                <!DOCTYPE html>
                <html>
                <head>
                    <style>
                    .row {
                        display : flex;
                        align-items : center;
                        margin-bottom: 15px;
                    }
                    .box {
                    height: 20px;
                    width: 20px;
                    border: 1px solid white;
                    margin-right : 5px;
                    }
                    .red {
                    background-color: #ef5777;
                    }
                    .blue {
                    background-color: #4bcffa;
                    }
                    .green {
                    background-color: #ffc048;
                    }
                    </style>
                </head>
                <body>
                    <div class="row">
                        <div class='box red'></div>
                        <span style="color:#ecf0f1"> $Topic1</span>
                    </div>
                    <div class="row">
                        <div class='box blue'></div>
                        <span style="color:#ecf0f1"> $Topic2</span>
                    </div>
                    <div class="row">
                        <div class='box green'></div>
                        <span style="color:#ecf0f1"> $Topic3</span>
                    </div>
                </body>
                </html>
                """
            )
        
        text = t.substitute({
            'Topic1': query_one, 
            'Topic2': query_two,
            'Topic3': 'Common MeSH Terms'})

        components.html(text, height=120)

        st.info('NOTE: Interact With Mouse.')
        st.warning('Use - Desktop Mode if on Mobile Devices')
        
        MKGraph.write_html("mesh_knowledge_graph.html")

        HtmlFile = open("mesh_knowledge_graph.html", 'r', encoding='utf-8')
        source_code = HtmlFile.read()
        components.html(source_code, height=1200, width=1650)    

try:
    os.remove("mesh_knowledge_graph.html")
except OSError:
    pass	

st.stop()     
