#!/usr/bin/env python
# coding: utf-8

# Fetching PubMed article metadata
from Bio import Entrez
from Bio import Medline

# Graph creation and visualisation
from pyvis.network import Network

def fetch_pmids(query_term, retrive='5'):
    """Search query term for corresponding pubmed article's pmids"""
    
    Entrez.email = 'akishirsath@gmail.com'
    
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax=retrive,
                            retmode='xml', 
                            term=query_term)
    
    results = Entrez.read(handle)
    
    return results

def fetch_articles(pmids, min_date='2016/01', max_date='2021/01'):
    """Returns pubmed data associated with pmids"""
    
    Entrez.email = 'akishirsath@gmail.com'

    handle = Entrez.efetch(db="pubmed", 
                           mindate=min_date, 
                           maxdate=max_date, 
                           id=pmids, 
                           rettype="medline", 
                           retmode="text")

    records = Medline.parse(handle)    
    
    return list(records)

def common_mesh_terms(articles_data_one, articles_data_two):    
    unique_terms_first = set()

    for article in articles_data_one:
        mesh_terms = article.get('MH', 'NO_MESH') 

        if mesh_terms != 'NO_MESH':
            for terms in mesh_terms:
                temp_list = [term.replace('*', '').replace(',', '').strip() for term in terms.split('/')]
                unique_terms_first.update(set(temp_list))

    unique_terms_sec = set()

    for article in articles_data_two:
        mesh_terms = article.get('MH', 'NO_MESH') 

        if mesh_terms != 'NO_MESH':
            for terms in mesh_terms:
                temp_list = [term.replace('*', '').replace(',', '').strip() for term in terms.split('/')]
                unique_terms_sec.update(set(temp_list))

    common_terms = (unique_terms_first & unique_terms_sec)

    return common_terms

def create_graph(articles_data, MKGraph, common_terms, color_type='FIRST'):
    
    if color_type == 'FIRST':
        colors = ["#ef5777", '#ffc048']
    if color_type == 'SECOND':
        colors = ["#4bcffa",'#ffc048']
    
    for article in articles_data:

        mesh_terms = article.get('MH', 'NO_MESH')

        if mesh_terms != 'NO_MESH':

            main_node = f"PMID_{article.get('PMID')}"
            MKGraph.add_node(main_node, size=35, title=article.get('TI'), color=colors[0])

            for terms in mesh_terms:
                temp_list = [term.replace('*', '').replace(',', '').strip() for term in terms.split('/')]

                primary_node = temp_list[0]
                secondary_nodes =  temp_list[1:]
                
                if primary_node not in common_terms:
                    MKGraph.add_node(primary_node, size=25, color=colors[0])
                else:
                    MKGraph.add_node(primary_node, size=25, color=colors[1])
                
                MKGraph.add_edge(main_node, primary_node)

                if len(secondary_nodes)>1:
                    for node in secondary_nodes:
                        if node not in common_terms:
                            MKGraph.add_node(node, size=15, color=colors[0])
                        else:
                            MKGraph.add_node(node, size=15, color=colors[1])
                        MKGraph.add_edge(primary_node, node)
        else:
            main_node = f"PMID_{article.get('PMID')}"
            primary_node = 'NO_MESH'

            MKGraph.add_node(main_node, size=35, title=article.get('TI'), color='#c56cf0')
            MKGraph.add_node(primary_node, size=25, title='NO MESH TERMS AVAILABLE', color='#c56cf0')

            MKGraph.add_edge(main_node, primary_node)    
    
    return MKGraph
