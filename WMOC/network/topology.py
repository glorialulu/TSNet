"""
The WMOC.network.topology figure out the tolopogy, i.e., 
upstream and downstream adjacent links for each pipe, and 
store the information in lists.

"""

def topology(wn, npipe):
    """Figure out the topology of the network
    
    Parameters
    ----------
    wn : wntr.network.model.WaterNetworkModel
        .inp file used for EPAnet simulation
    npipe : integer
        Number of pipes

    Returns
    -------
    links1 : list 
        The id of ajacent pipe on the start node. 
        The sign represents the direction of the pipe. 
        + : flowing into the junction
        - : flowing out from the junction
    links2 : list 
        The id of ajacent pipe on the end node. 
        The sign represents the direction of the pipe. 
        + : flowing into the junction
        - : flowing out from the junction
    utype : list 
        The type of the upstream ajacent links. 
        If the link is not pipe, the name of that link 
        will also be included.
        If there is no upstream link, the type of the start node 
        will be recorded.  
    dtype : list
        The type of the downstream ajacent links. 
        If the link is not pipe, the name of that link 
        will also be included.
        If there is no downstream link, the type of the end node 
        will be recorded.
    """
    G = wn.get_graph()
    
    length = wn.query_link_attribute('length')
    G.weight_graph(link_attribute = length)

    # add 'id' attribute to networkx links 
    i =1
    for ln, link in wn.links():          
        G.edges[link.start_node_name,link.end_node_name, ln]['id'] = i
        i+=1 
        
    # allocate the parameters
    links1 = [0] * len(wn.links) 
    links2 = [0] * len(wn.links)  
    utype = [('Pipe',0)] * npipe
    dtype = [('Pipe',0)] * npipe
 
    # Adjcant pipes for each pipe IN:+; OUT:-   
    for _, link in wn.links():
        pn = link.id
        links1[int(pn)-1] = [int(p['id'])
                           for _, attr in G.pred[link.start_node_name].items() 
                           for _,p in attr.items()
                           if p['id'] != pn] 
        
        for _, attr in G.succ[link.start_node_name].items() :
            for _,p in attr.items():
                if p['id'] != pn:
                    links1[int(pn)-1].append(-1* int(p['id']))
               
                
        # right (end) adjcant pipes 
        links2[int(pn)-1] = [int(p['id']) 
                       for _, attr in G.pred[link.end_node_name].items() 
                       for _,p in attr.items()
                       if p['id'] != pn]  
        
        for _, attr in G.succ[link.end_node_name].items() :
            for _,p in attr.items():
                if p['id'] != pn:
                    links2[int(pn)-1].append(-1*int(p['id']))
                    
    #figure out downstream type and upstream type                
    for _,pipe in wn.pipes():
        pn = pipe.id-1 
        if links1[pn] :   
            if max(map(abs, links1[pn])) > npipe:   
                utype[pn] = [(l.link_type,l.name)
                         for _,l in wn.links() 
                         if l.id == abs(links1[pn][0])][0] 
                
                if links1[abs(links1[pn][0])-1] and links2[abs(links1[pn][0])-1]:
                    links1[pn] = [i 
                    for i in [links1[abs(links1[pn][0])-1], links2[abs(links1[pn][0])-1]]
                    if abs(i[0]) -1 != pn][0]         
                else:
                    links1[pn] = ['End']
                    
        else:
            utype[pn] = (wn.nodes[pipe.start_node_name].node_type, 
                         wn.nodes[pipe.start_node_name])

        if links2[pn] :            
            if max(map(abs, links2[pn])) > npipe:
                dtype[pn] = [(l.link_type,l.name) 
                         for _,l in wn.links() 
                         if l.id == abs(links2[pn][0])][0]   
                 
                if links1[abs(links2[pn][0])-1] and links2[abs(links2[pn][0])-1]:
                    links2[pn] = [i 
                    for i in [links1[abs(links2[pn][0])-1], links2[abs(links2[pn][0])-1]]
                    if abs(i[0]) -1 != pn][0]  
                else:
                    links2[pn] = ['End']
                           
        else:
            dtype[pn] = (wn.nodes[pipe.end_node_name].node_type,
                         wn.nodes[pipe.end_node_name])
        
    return links1, links2, utype, dtype