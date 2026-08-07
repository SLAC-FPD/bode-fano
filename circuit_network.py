from circuit_reader import *
from circuit_calcs import *

cd = CircuitData() # create instance.
# options
template = "bfe"
param_file = "params/bfe_params_251218.txt"
variation = "parallel_bfe"

# initialize
cd.read_template(template, param_file, variation)

def make_node_dict(component_dict):
    # read circuit_text and get components connected to nodes
    node_dict = {}  # let's make this to begin with
    # find nodes corresponding to where
    for component_key in component_dict.keys():
        component = component_dict[component_key]
        for node in component.nodes:
            try: node_dict[node].append(component_key)  # add to node
            except KeyError: node_dict[node] = [component_key]  # new node!
    # print(node_dict)
    return node_dict

def connect_components(component_dict, node_dict, component_name, start_node=None, ground_node='0'):
    new_component_dict = component_dict.copy()
    new_node_dict = node_dict.copy()
    component = new_component_dict[component_name]  # start component
    if start_node is None:
        if ground_node in component.nodes: start_node = ground_node
        else: start_node = component.nodes[0]  # whatever first component we have
    connected_node = component.nodes[1] if component.nodes[0] is str(start_node) else str(component.nodes[0])  # find connected node
    del new_component_dict[component.label]  # remove that component in new_component_dict...
    new_node_dict[start_node].remove(component.label)  # and new_node_dict.
    new_node_dict[connected_node].remove(component.label)  # and new_node_dict (again)
    connected_components = new_node_dict[connected_node]
    print(component_name, start_node, connected_node)
    while len(connected_components) > 0 and connected_node != ground_node:
        connected_component_name = connected_components[0]
        start_node = connected_node
        if connected_node == ground_node: continue
        connected_component = new_component_dict[connected_component_name]
        component.comps_to[connected_component.label] = None
        connected_component.comps_from[component.label] = None
        connect_components(new_component_dict, new_node_dict, connected_component.label, start_node, ground_node)  # keep connecting
    # while connected_node != ground_node and start_node != connected_node:
    #     start_node = connected_node
    #     print("NODES", start_node, connected_node, ground_node)
    #     print("CON COMPS:", connected_components)
    #     print(new_component_dict)
    #     print(new_node_dict)

def make_network(component_dict, component_name, start_node=None, ground_node='0'):
    node_dict = make_node_dict(component_dict)
    connect_components(component_dict, node_dict, component_name, start_node, ground_node)

def reset_intermediate_impedances(component_dict):
    print(component_dict.keys())
    for component_key in component_dict.keys():
        component = component_dict[component_key]
        component.reset_intermediate_impedance()

def t_circuitify(component_dict):
    new_component_dict = component_dict.copy()
    t_node = 0
    node_changes = {}
    # t_circuitifies all mutual inductance in a component dictionary. In other words, removes all k's and replaces with new l1, l2, m's
    for component_key in component_dict.keys():  # iterate over old dictionary, add to new dictionary.
        if component_key[0] == "k":  # Let's make that t-circuit!
            k = component_dict[component_key]  # get the mutual component
            del new_component_dict[component_key]  # remove it from new dictionary
            l0 = component_dict[k.nodes[0]]
            l1 = component_dict[k.nodes[1]]
            if l0.nodes[0] != '0' : l0_node_0, l0_node_1 = l0.nodes[0], l0.nodes[1]
            else: l0_node_0, l0_node_1 = l0.nodes[1], l0.nodes[0]
            if l1.nodes[0] != '0' : l1_node_0, l1_node_1 = l1.nodes[0], l1.nodes[1]
            else: l1_node_0, l1_node_1 = l1.nodes[1], l1.nodes[0]
            if l0.nodes[1] != '0': node_to_keep, node_to_change = l0_node_1, l1_node_1  # as long as it's not the ground node
            else: node_to_keep, node_to_change = l1_node_1, l0_node_1
            node_to_add = f"t{t_node}"  # the t-type node
            m = Component()
            m.value = k.value*np.sqrt(l0.value*l1.value)
            m.label = f"l{l0.label[1:]}_{l1.label[1:]}"  # combined two names
            if m.label in component_dict.keys(): m.label += "_"  # if label name is redundant...
            m.nodes = [node_to_add, node_to_keep]
            l0.value = l0.value - m.value
            l0.nodes = [l0_node_0, node_to_add]
            l1.value = l1.value - m.value
            l1.nodes = [node_to_add, l1_node_0]
            new_component_dict[m.label] = m  # add to new dictionary
            node_changes[node_to_change] = node_to_keep
            print(l0_node_0, l0_node_1, l1_node_0, l1_node_1)
            t_node += 1
        else: continue
    for node_to_change in node_changes.keys():
        node_to_keep = node_changes[node_to_change]
        for component_key in new_component_dict.keys():  # combine nodes...
            if node_to_change in new_component_dict[component_key].nodes:
                new_component_dict[component_key].nodes = [node_to_keep if x == node_to_change else x for x in new_component_dict[component_key].nodes]
    return new_component_dict

class Component:  # class that separates each component in a .cir file
    def __init__(self, circuit_text=None):
        self.reset()
        if circuit_text is not None: self.read_component(circuit_text)
        # print(self)
        
    def reset(self):
        self.label = None  # label for values e.g., leff
        self.value_label = None
        self.value = None  # magnitude of component (resistance, inductance, capacitance, coupling, critical current)
        self.nodes = None  # nodes of a .cir file
        self.type = None  # r=resistance, l=inductance, c=capacitance, k=mutual, b=junction
        self.other = None  # other miscellaneous things that might be relevant later
        self.supported = True  # is this a supported component?
        # below are network specific things and are directional!
        self.impedance = None
        self.intermediate_impedance = None  # the directional intermediate impedance thing
        self.comps_from = {}  # if one, in series, if multiple, in parallel.
        self.comps_to = {}
    
    def __str__(self):
        return f"Component: {self.label}\nNodes: {self.nodes}\nValue Label: {self.value_label}\nValue: {self.value}"
    
    def get_impedance(self, freq, phase=0):
        # gives impedance of current branch and sets it as self impedance.
        # only use phase if jj, can be fed separately
        omega = freq * 2 * np.pi
        if self.type is None:
            try: self.type = self.label[0]
            except: impedance = 0
        if self.type == "r":
            impedance = self.value
        elif self.type == "l":
            impedance = self.value * omega * 1j
        elif self.type == "c":
            impedance = 1 / (self.value * omega * 1j)
        elif self.type == "b":  # maybe turn into complex later
            lj = calc_lj(self.value, phase)
            impedance = lj * omega * 1j
        else: impedance = 0
        self.impedance = impedance
        return impedance
    
    def reset_intermediate_impedance(self):  # just this one
        for comps_key in self.comps_from.keys(): self.comps_from[comps_key] = None
        for comps_key in self.comps_to.keys(): self.comps_to[comps_key] = None
    
    def read_component(self, circuit_text):  # read from circuit text line
        split_text = circuit_text.strip().split(" ")
        self.label = split_text[0]
        measurable_type = self.label[0]  # first letter of label
        self.type = measurable_type
        if self.type in "rlc":  # current extent of support
            # note: r may need to consider m=1, etc.
            self.nodes = split_text[1:3]  # first two values after name
            self.value_label = split_text[3][1:-1]  # remove curly braces...
            try: self.other = split_text[4:]  # dump everything else, if it exists
            except: pass
        elif self.type == "k":
            self.nodes = split_text[1:3]
            self.value_label = split_text[3][1:-1]  # remove curly braces...
        elif self.type in "biv":  # i, v may need some more adjusting
            self.nodes = []  # define the empty list
            self.other = []  # define the empty list
            for entry_text in split_text:
                try:
                    int(entry_text)
                    self.nodes.append(entry_text)  # if it is an int
                except:
                    if 'ics' in entry_text or f'{self.label}_mag' in entry_text:  # critical current assigned to value_label
                        self.value_label = entry_text.split("=")[-1][1:-1]  # lol
                    else: self.other.append(entry_text)
        else: self.supported = False
        try: self.value = float(self.value_label)
        except ValueError or TypeError: pass


def get_intermediate_impedances(component_dict, freq, phase, end_components=None, ground_node='0'):
    if end_components is None:  # find end components. intermediate impedance is self impedance when first updated
        end_components = []
        for component_key in component_dict.keys():
            component = component_dict[component_key]
            if len(component.comps_to.keys()) < 1:
                end_components.append(component.label)
                component.intermediate_impedance = component.get_impedance(freq, phase)  # dead end gets own self impedance
           
    new_end_components = end_components.copy()
    for end_component_name in end_components:
        print("END COMPS: ", end_components, new_end_components)
        print("CUR COMP: ", end_component_name)
        # currently assumes only one from component. maybe need to see how to do multiple later on.
        end_component = component_dict[end_component_name]
        if end_component.intermediate_impedance is None: end_component.intermediate_impedance = end_component.get_impedance(freq, phase)
        if end_component_name not in new_end_components: break  # removed via parallel
        from_components = list(end_component.comps_from.keys())  # need to turn into list first
        print("FROM COMPS: ", from_components, len(from_components))
        if len(from_components) > 0:  # we're not at the start
            from_component_name = from_components[0]  # support multiple from components later?
            from_component = component_dict[from_component_name]
            from_to_components = list(from_component.comps_to.keys())  # need to turn into list first
            if len(from_to_components) == 1:  # it's in series
                from_component.intermediate_impedance = end_component.intermediate_impedance + from_component.get_impedance(freq, phase)
                new_end_components.remove(end_component_name)
                new_end_components.append(from_component_name)  # update end component
            elif len(from_to_components) >= 1:  # it's in parallel
                new_from_to_components = from_to_components.copy()
                for from_to_component_name in from_to_components:
                    if from_to_component_name in end_components: new_from_to_components.remove(from_to_component_name)
                if len(new_from_to_components) == 0:  # only if ALL from_to_components have been removed, we can do a parallel calculations.
                    intermediate_impedance_inverted = 0  # do parallel calcs
                    for from_to_component_name in from_to_components:  # reiterate...
                        from_to_component = component_dict[from_to_component_name]
                        print("PARALLELS", from_to_component.label, from_to_component.impedance)
                        intermediate_impedance_inverted += 1 / from_to_component.intermediate_impedance
                        new_end_components.remove(from_to_component_name)
                    from_component.intermediate_impedance = from_component.get_impedance(freq, phase) + 1 / intermediate_impedance_inverted  # inverse after adding all parallel branches
                    new_end_components.append(from_component_name)  # update end component
                else: continue
            print("INT IMP CALC RESULT: ", from_component_name, from_component.intermediate_impedance)
            get_intermediate_impedances(component_dict, freq, phase, new_end_components)
        else: break  # no more from components, we're at the start!


comp_dict = {}
for comp_key in cd.circuit_text.keys():
    comp_key_text = cd.circuit_text[comp_key]
    if comp_key[0] == "k" and comp_key_text[0] == "l": # we have a correct? mutual inductance
        pass  # "passes" the test
    elif comp_key[0] in "rlcbiv":
        try: int(comp_key_text[0]) # see if we have a correct component that has a "node"
        except ValueError: continue
    else: continue
    circuit_text = f"{comp_key} {comp_key_text}"
    comp_dict[comp_key] = Component(circuit_text)
    comp_dict[comp_key].value = cd.params[comp_dict[comp_key].value_label]  # this seems awfully roundabout
    # print(comp_dict[comp_key])
# print(comp_dict)
comp_dict = t_circuitify(comp_dict)

make_network(comp_dict, "v1")
for comp_key in comp_dict.keys():
    comp = comp_dict[comp_key]
    print(comp)
    print("FROM", comp.comps_from)
    print("TO", comp.comps_to, "\n")
get_intermediate_impedances(comp_dict, 5e7, 0, ['l3', 'rout', 'l1_2'])  # 'rin'
