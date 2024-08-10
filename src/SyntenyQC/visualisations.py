# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 16:31:23 2024

@author: u03132tk
"""
import networkx as nx
import plotly.graph_objects as go
import logging

def write_graph(graph : nx.Graph, path : str, similarity_filter : float, 
                min_edge_view : float, logger_name : str) -> None:
    '''
    Write graph to html figure.

    Parameters
    ----------
    graph : nx.Graph
        Graph to write.
    path : str
        Path to write graph to.
    similarity_filter : float
        Minimum similarity for two nodes to be considered equivalent - used for 
        edge colouring.
    min_edge_view : float
        Minimum similarity for an edge between two nodes to be shown - used for 
        edge colouring.
    logger_name : str
        name of logger to record results
    '''
    #define node positions 
    pos = define_pos(graph)
    nx.set_node_attributes(graph, pos, 'pos')
    
    #build traces
    edge_traces = make_edge_traces(graph, similarity_filter, min_edge_view)
    node_trace = make_node_traces(graph)
    
    #build figure
    title = f'Reciprocal Best Hits network (similarity filter {similarity_filter}, '\
                 'min edge shown {min_edge_view})'
    fig = go.Figure(data=[*edge_traces, 
                          node_trace],
                    layout=go.Layout(title=title,
                                     titlefont_size=16,
                                     showlegend=False,
                                     hovermode='closest',
                                     margin=dict(b=20,l=5,r=5,t=40),
                                     xaxis=dict(showgrid=False, 
                                                zeroline=False, 
                                                showticklabels=False),
                                     yaxis=dict(showgrid=False, 
                                                zeroline=False, 
                                                showticklabels=False),
                                     paper_bgcolor='rgba(0,0,0,0)',
                                     plot_bgcolor='rgba(0,0,0,0)'
                                     )
                    )
    
    #write fig to html file
    fig.write_html(path)
    
    #log details
    logger = logging.getLogger(logger_name)
    logger.info(f'Made RBH graph of {len(graph.nodes)} unpruned '\
                    f'neighbourhoods - written to {path}')
    
    
def define_pos(graph : nx.Graph) -> dict:
    '''
    Define node positions using spring layout.  Note, for visual clarity position 
    is defined based on a copy of the input graph where 
    new_edge = (graph_edge ** 2) * 1.5.  

    Parameters
    ----------
    graph : nx.Graph
        Graph to get pos from.

    Returns
    -------
    pos : dict
        Position mapping for each node.

    '''    
    all_nodes = graph.nodes
    copy_graph = nx.Graph()
    copy_graph.add_nodes_from(all_nodes)
    edges = [(n1, n2, (e**2)*1.5) for n1, n2, e in graph.edges.data('weight')]
    copy_graph.add_weighted_edges_from(edges)
    pos = nx.drawing.layout.spring_layout(graph)
    return pos
    
def make_edge_traces(graph : nx.Graph, similarity_filter : float, 
                     min_edge_view : float) -> list:
    '''
    Make edge traces. 

    Parameters
    ----------
    graph : nx.Graph
        Graph to build traces from.
    similarity_filter : float
        Minimum similarity for two nodes to be considered equivalent - used for 
        edge colouring.
    min_edge_view : float
        Minimum similarity for an edge between two nodes to be shown - used for 
        edge colouring.

    Returns
    -------
    edges : list
        List of edge traces.

    '''
    
    edges = []
    
    #for every edges
    for n1, n2, weight in graph.edges.data("weight"):
        
        #define edge colour 
        if weight >= similarity_filter and min_edge_view < similarity_filter:
            #red if it is above threshold
            colour = f'rgba(255,0,0,{weight})'
        else:
            #black if under threshold
            colour = f'rgba(0,0,0,{weight})'
        
        #def positions of edge termini
        x0, y0 = graph.nodes[n1]['pos']
        x1, y1 = graph.nodes[n2]['pos']
        x=[x0, x1, None]
        y=[y0,y1,None]
        
        #make edge trace for line
        edges += [go.Scatter(x=x, 
                             y=y,
                             line=dict(width = (weight**2)*1.5, 
                                       color = colour
                                       ),
                             hoverinfo='none',
                             mode='lines'
                             )
                  ]
                  
        #make edge trace for the edge text
        middle_x = [(x0+x1)/2]
        middle_y = [(y0+y1)/2]
        middle_text = [f'{n1} - {n2} - {round(weight,2)}']
        edges += [go.Scatter(x=middle_x,
                             y=middle_y,
                             text=middle_text,
                             mode='markers',
                             hoverinfo='text',
                             marker=go.Marker(opacity=0
                                              )
                             )
                  ]
    return edges


def make_node_traces(graph : nx.Graph) -> go.Scatter:
    '''
    Make node traces

    Parameters
    ----------
    graph : nx.Graph
        Graph with nodes of interest.

    Returns
    -------
    node_trace : go._scatter.Scatter
        A node trace.

    '''
    #get node position 
    node_x = []
    node_y = []
    for node in graph.nodes():
        x, y = graph.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)
    
    #build node trace 
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))
    
    #build node annotations
    node_degrees = []
    for adjacencies in graph.adjacency():
        node_degrees.append(len(adjacencies[1]))
    node_text = []
    for degree, node in zip(node_degrees, graph.nodes):
        node_text.append(f"{node} # of connections: {degree}")
    
    #colour nodes by degree and annotate each node
    node_trace.marker.color = node_degrees
    node_trace.text = node_text
    return node_trace


def write_hist(graph : nx.Graph, path : str, logger_name : str) -> None:
    '''
    Write histogram for edge weights in graph.

    Parameters
    ----------
    G : nx.Graph
        Graph with edges of interest.
    path : str
        Path to write histogram.
    logger_name : str
        name of logger to write results
    '''
    
    #get edge weights
    all_weights = [weight for n1, n2, weight in graph.edges.data("weight")]
    
    #make hist figure and write to html
    fig = go.Figure(data=[go.Histogram(x=all_weights,
                                       
                                       xbins=dict( # bins used for histogram
                                                    start=0,
                                                    end=1,
                                                    size=0.05
                                                ),
                                       autobinx = False
                                                )],
                    layout=go.Layout(
                           title='RBH similarities',
                           titlefont_size=16,
                           showlegend=False,
                           hovermode='closest',
                           margin=dict(b=20,l=5,r=5,t=40),
                           xaxis=dict (range= [0, 1]),
                           paper_bgcolor='rgba(0,0,0,0)',
                           plot_bgcolor='rgba(0,0,0,0)'
                           ))
    fig.write_html(path)
    
    #log results
    logger = logging.getLogger(logger_name)
    logger.info('Made histogram showing the distribution of RBH similarities '\
                    f'between all {len(graph.nodes)} neighbourhoods - written '\
                        f'to {path}')