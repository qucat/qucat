from circuit_elements import L,J,C,R,Series,Parallel,Component

class Node(object):
    def __init__(self, label):
        self.parent = None
        self.label = label
    def is_equal_to(self,other_node):
        return self.label == other_node.label
class ConnectionNode(Node):
    def __init__(self, label,child_left,child_right):
        super(ConnectionNode,self).__init__(label)
        self.child_left = child_left
        self.child_right = child_right
        self.set_child_parenthood()
        
    def set_child_parenthood(self):
        self.child_left.parent = self
        self.child_right.parent = self

class SeriesNode(ConnectionNode):
    """docstring for SeriesNode"""
    def __init__(self, label,child_left,child_right):
        super(SeriesNode, self).__init__(label,child_left,child_right)

class ParallelNode(ConnectionNode):
    """docstring for ParallelNode"""
    def __init__(self, label,child_left,child_right):
        super(ParallelNode, self).__init__(label,child_left,child_right)
        
class ComponentNode(Node):
    def __init__(self, label,component):
        super(ComponentNode,self).__init__(label)
        self.component = component
        
            
def circuit_to_tree(circuit):
    global label
    label = 0
    component_node_dict = {}
    
    def recursive_function(circuit):
        global label
        label +=1
        if type(circuit) == Parallel:
            return ParallelNode(
                label = label,
                child_left = recursive_function(circuit.component_left),
                child_right = recursive_function(circuit.component_right))
        elif type(circuit) == Series:
            return SeriesNode(
                label = label,
                child_left = recursive_function(circuit.component_left),
                child_right = recursive_function(circuit.component_right))
        elif isinstance(circuit,Component):
            cn = ComponentNode(
                label = label,
                component = circuit)
            component_node_dict[circuit.label] = cn
            return cn
            
    return recursive_function(circuit),component_node_dict

def tree_to_circuit(tree):
    if type(tree) == SeriesNode:
        return Series(tree_to_circuit(tree.child_left),
                   tree_to_circuit(tree.child_right))
    elif type(tree) == ParallelNode:
        return Parallel(tree_to_circuit(tree.child_left),
                   tree_to_circuit(tree.child_right))
    elif type(tree) == ComponentNode:
        return tree.component


def rotate_tree(rotation_node):
    
    # Check if the tree is already correctly rotated
    if rotation_node.parent.parent == None:
        return rotation_node.parent
    
    def recursive_function(node):
        parent_node = node.parent
        if parent_node.parent == None:
            if parent_node.child_left.is_equal_to(node):
                return parent_node.child_right
            elif parent_node.child_right.is_equal_to(node):
                return parent_node.child_left
        else:
            if parent_node.child_left.is_equal_to(node):
                parent_node.child_left = recursive_function(parent_node)
            elif parent_node.child_right.is_equal_to(node):
                parent_node.child_right = recursive_function(parent_node)
            return parent_node
        
    rotated_tree = ParallelNode(label=1,
                        child_left = rotation_node,
                        child_right = recursive_function(rotation_node))
    
    # Reset all node parenthoods
    def set_parenthoods(node):
        if type(node) == ConnectionNode:
            node.set_child_parenthood()
            set_parenthoods(node.child_left)
            set_parenthoods(node.child_right)
    
    set_parenthoods(rotated_tree)
    return rotated_tree

def rotate_circuit(circuit,rotation_component):
    tree,components = circuit_to_tree(circuit)
    rotated_tree = rotate_tree(components[rotation_component])
    return tree_to_circuit(rotated_tree)