#!/usr/bin/env python3

from ete3 import Tree
import pandas as pd

def parse_tree_topology_carefully():
    """
    Carefully parse the phylogenetic tree to understand true topology
    """
    print("=== Detailed Tree Topology Analysis ===")
    
    # Load tree
    tree_file = '/scratch1/chris/introner_vis/iqtree/mpusilla.snps.4d.notMT.min4.phy.treefile'
    tree = Tree(tree_file)
    
    # Root with outgroups
    outgroups = ['RCC1749', 'RCC3052']
    outgroup_nodes = []
    for node in tree.traverse():
        if node.is_leaf() and node.name in outgroups:
            outgroup_nodes.append(node)
    
    if len(outgroup_nodes) == 2:
        outgroup_ancestor = tree.get_common_ancestor(outgroup_nodes)
        tree.set_outgroup(outgroup_ancestor)
        print("Tree rooted with outgroups RCC1749, RCC3052")
    
    # Define introner groups
    introner_group1 = ['RCC1614', 'RCC1698', 'RCC373', 'RCC465']
    introner_group2 = ['CCMP1545', 'RCC114', 'RCC2482', 'RCC629', 'RCC692', 'RCC693', 'RCC833']
    
    print("\nIntroner Group Definitions:")
    print(f"Group1: {introner_group1}")
    print(f"Group2: {introner_group2}")
    
    # Analyze actual tree structure
    print("\n=== DETAILED TREE STRUCTURE ANALYSIS ===")
    
    # Get all Group1 samples and find their relationships
    group1_nodes = {}
    group2_nodes = {}
    
    for node in tree.traverse():
        if node.is_leaf():
            if node.name in introner_group1:
                group1_nodes[node.name] = node
            elif node.name in introner_group2:
                group2_nodes[node.name] = node
    
    print(f"\nFound Group1 samples in tree: {list(group1_nodes.keys())}")
    print(f"Found Group2 samples in tree: {list(group2_nodes.keys())}")
    
    # Test monophyly more carefully
    print("\n=== MONOPHYLY TESTING ===")
    
    # Test Group1 monophyly
    if len(group1_nodes) > 1:
        group1_list = list(group1_nodes.values())
        group1_ancestor = tree.get_common_ancestor(group1_list)
        group1_descendants = [leaf.name for leaf in group1_ancestor.get_leaves()]
        
        # Check if only Group1 samples are descendants
        non_group1_in_clade = [name for name in group1_descendants if name not in introner_group1 and name not in outgroups]
        
        print(f"\nGroup1 Analysis:")
        print(f"  Common ancestor descendants: {group1_descendants}")
        print(f"  Non-Group1 samples in clade: {non_group1_in_clade}")
        print(f"  Group1 is monophyletic: {len(non_group1_in_clade) == 0}")
    
    # Test Group2 monophyly
    if len(group2_nodes) > 1:
        group2_list = list(group2_nodes.values())
        group2_ancestor = tree.get_common_ancestor(group2_list)
        group2_descendants = [leaf.name for leaf in group2_ancestor.get_leaves()]
        
        # Check if only Group2 samples are descendants
        non_group2_in_clade = [name for name in group2_descendants if name not in introner_group2 and name not in outgroups]
        
        print(f"\nGroup2 Analysis:")
        print(f"  Common ancestor descendants: {group2_descendants}")
        print(f"  Non-Group2 samples in clade: {non_group2_in_clade}")
        print(f"  Group2 is monophyletic: {len(non_group2_in_clade) == 0}")
    
    # Let's examine the relationships more systematically
    print("\n=== SYSTEMATIC CLADE ANALYSIS ===")
    
    # Find the smallest clades containing specific samples
    target_samples = introner_group1 + introner_group2
    
    for sample in target_samples:
        if sample in [node.name for node in tree.traverse() if node.is_leaf()]:
            sample_node = None
            for node in tree.traverse():
                if node.is_leaf() and node.name == sample:
                    sample_node = node
                    break
            
            if sample_node:
                # Find sister groups
                parent = sample_node.up
                if parent:
                    sisters = [child for child in parent.get_children() if child != sample_node]
                    sister_names = []
                    for sister in sisters:
                        if sister.is_leaf():
                            sister_names.append(sister.name)
                        else:
                            sister_names.extend([leaf.name for leaf in sister.get_leaves()])
                    
                    print(f"{sample}: sisters = {sister_names}")
    
    # Let's also check specific pairwise relationships mentioned
    print("\n=== SPECIFIC RELATIONSHIPS ===")
    
    samples_to_check = ['CCMP1545', 'RCC692', 'RCC114', 'RCC833', 'RCC2482']
    
    for i, sample1 in enumerate(samples_to_check):
        for sample2 in samples_to_check[i+1:]:
            # Find nodes
            node1 = None
            node2 = None
            for node in tree.traverse():
                if node.is_leaf():
                    if node.name == sample1:
                        node1 = node
                    elif node.name == sample2:
                        node2 = node
            
            if node1 and node2:
                # Find common ancestor
                common_ancestor = tree.get_common_ancestor([node1, node2])
                descendants = [leaf.name for leaf in common_ancestor.get_leaves()]
                other_descendants = [name for name in descendants if name not in [sample1, sample2]]
                
                print(f"{sample1} + {sample2}: common ancestor includes {other_descendants}")
    
    # Print tree in ASCII for visual inspection
    print("\n=== TREE TOPOLOGY (ASCII) ===")
    try:
        print(tree.get_ascii(show_internal=True))
    except:
        print("ASCII tree display not available")
    
    # Return detailed topology analysis
    return tree, group1_nodes, group2_nodes

def reanalyze_introner_groups_based_on_tree(tree, group1_nodes, group2_nodes):
    """
    Reanalyze introner groups based on actual tree topology
    """
    print("\n=== REVISED POPULATION STRUCTURE ANALYSIS ===")
    
    # Let's identify true monophyletic groups in the tree
    all_samples = list(group1_nodes.keys()) + list(group2_nodes.keys())
    
    # Try to identify natural monophyletic groups
    print("Identifying natural monophyletic clades...")
    
    # Get all internal nodes and check what they contain
    clades_found = []
    
    for node in tree.traverse():
        if not node.is_leaf() and len(node.get_leaves()) > 1:
            descendants = [leaf.name for leaf in node.get_leaves()]
            # Filter to only our samples of interest
            relevant_descendants = [name for name in descendants if name in all_samples]
            
            if len(relevant_descendants) > 1:
                # Check if this is a "clean" clade (doesn't contain outgroups or obsolete samples)
                clean_clade = True
                for desc in descendants:
                    if desc in ['RCC1749', 'RCC3052', 'CCMP490', 'RCC835', 'RCC647']:
                        clean_clade = False
                        break
                
                if clean_clade and len(relevant_descendants) >= 2:
                    clades_found.append({
                        'size': len(relevant_descendants),
                        'samples': relevant_descendants,
                        'all_descendants': descendants
                    })
    
    # Sort clades by size
    clades_found.sort(key=lambda x: x['size'])
    
    print(f"\nMonophyletic clades found:")
    for i, clade in enumerate(clades_found):
        print(f"  Clade {i+1} (n={clade['size']}): {clade['samples']}")
    
    # Compare with introner groupings
    introner_group1 = ['RCC1614', 'RCC1698', 'RCC373', 'RCC465']
    introner_group2 = ['CCMP1545', 'RCC114', 'RCC2482', 'RCC629', 'RCC692', 'RCC693', 'RCC833']
    
    print(f"\nComparison with introner groups:")
    
    for i, clade in enumerate(clades_found):
        group1_overlap = len(set(clade['samples']) & set(introner_group1))
        group2_overlap = len(set(clade['samples']) & set(introner_group2))
        
        print(f"  Clade {i+1}: {group1_overlap}/4 Group1, {group2_overlap}/7 Group2")
        
        if group1_overlap == len(introner_group1) and group1_overlap == len(clade['samples']):
            print(f"    *** Clade {i+1} EXACTLY matches Introner Group1 ***")
        elif group2_overlap == len(introner_group2) and group2_overlap == len(clade['samples']):
            print(f"    *** Clade {i+1} EXACTLY matches Introner Group2 ***")

if __name__ == "__main__":
    tree, group1_nodes, group2_nodes = parse_tree_topology_carefully()
    reanalyze_introner_groups_based_on_tree(tree, group1_nodes, group2_nodes)