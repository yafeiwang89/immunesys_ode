 
# This file is auto-generated from a Python script that parses a PhysiCell configuration (.xml) file.
#
# Edit at your own risk.
#
import os
from ipywidgets import Label,Text,Checkbox,Button,HBox,VBox,FloatText,IntText,BoundedIntText,BoundedFloatText,Layout,Box
    
class UserTab(object):

    def __init__(self):
        
        micron_units = Label('micron')   # use "option m" (Mac, for micro symbol)

        constWidth = '180px'
        tab_height = '500px'
        stepsize = 10

        #style = {'description_width': '250px'}
        style = {'description_width': '25%'}
        layout = {'width': '400px'}

        name_button_layout={'width':'25%'}
        widget_layout = {'width': '15%'}
        units_button_layout ={'width':'15%'}
        desc_button_layout={'width':'45%'}

        param_name1 = Button(description='random_seed', disabled=True, layout=name_button_layout)
        param_name1.style.button_color = 'lightgreen'

        self.random_seed = IntText(
          value=0,
          step=1,
          style=style, layout=widget_layout)

        param_name2 = Button(description='infection_o2_relative_uptake', disabled=True, layout=name_button_layout)
        param_name2.style.button_color = 'tan'

        self.infection_o2_relative_uptake = FloatText(
          value=0.5,
          step=0.1,
          style=style, layout=widget_layout)

        param_name3 = Button(description='infection_persistence_time', disabled=True, layout=name_button_layout)
        param_name3.style.button_color = 'lightgreen'

        self.infection_persistence_time = FloatText(
          value=5.0,
          step=0.1,
          style=style, layout=widget_layout)

        param_name4 = Button(description='infection_migration_speed', disabled=True, layout=name_button_layout)
        param_name4.style.button_color = 'tan'

        self.infection_migration_speed = FloatText(
          value=0.4,
          step=0.1,
          style=style, layout=widget_layout)

        param_name5 = Button(description='unattached_infection_migration_bias', disabled=True, layout=name_button_layout)
        param_name5.style.button_color = 'lightgreen'

        self.unattached_infection_migration_bias = FloatText(
          value=0.0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name6 = Button(description='infection_apoptosis_rate', disabled=True, layout=name_button_layout)
        param_name6.style.button_color = 'tan'

        self.infection_apoptosis_rate = FloatText(
          value=0.0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name7 = Button(description='infection_relative_cycle_entry_rate', disabled=True, layout=name_button_layout)
        param_name7.style.button_color = 'lightgreen'

        self.infection_relative_cycle_entry_rate = FloatText(
          value=0.0166,
          step=0.001,
          style=style, layout=widget_layout)

        param_name8 = Button(description='scout_o2_relative_uptake', disabled=True, layout=name_button_layout)
        param_name8.style.button_color = 'tan'

        self.scout_o2_relative_uptake = FloatText(
          value=0.02,
          step=0.001,
          style=style, layout=widget_layout)

        param_name9 = Button(description='scout_apoptosis_rate', disabled=True, layout=name_button_layout)
        param_name9.style.button_color = 'lightgreen'

        self.scout_apoptosis_rate = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name10 = Button(description='scout_persistence_time', disabled=True, layout=name_button_layout)
        param_name10.style.button_color = 'tan'

        self.scout_persistence_time = FloatText(
          value=10,
          step=1,
          style=style, layout=widget_layout)

        param_name11 = Button(description='scout_migration_speed', disabled=True, layout=name_button_layout)
        param_name11.style.button_color = 'lightgreen'

        self.scout_migration_speed = FloatText(
          value=0.1,
          step=0.01,
          style=style, layout=widget_layout)

        param_name12 = Button(description='unattached_scout_migration_bias', disabled=True, layout=name_button_layout)
        param_name12.style.button_color = 'tan'

        self.unattached_scout_migration_bias = FloatText(
          value=0.7,
          step=0.1,
          style=style, layout=widget_layout)

        param_name13 = Button(description='scout_relative_cycle_entry_rate', disabled=True, layout=name_button_layout)
        param_name13.style.button_color = 'lightgreen'

        self.scout_relative_cycle_entry_rate = FloatText(
          value=0.0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name14 = Button(description='scout_relative_adhesion', disabled=True, layout=name_button_layout)
        param_name14.style.button_color = 'tan'

        self.scout_relative_adhesion = FloatText(
          value=0.0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name15 = Button(description='scout_relative_repulsion', disabled=True, layout=name_button_layout)
        param_name15.style.button_color = 'lightgreen'

        self.scout_relative_repulsion = FloatText(
          value=10.0,
          step=1,
          style=style, layout=widget_layout)

        param_name16 = Button(description='number_of_infection', disabled=True, layout=name_button_layout)
        param_name16.style.button_color = 'tan'

        self.number_of_infection = IntText(
          value=20,
          step=1,
          style=style, layout=widget_layout)

        param_name17 = Button(description='number_of_scouts', disabled=True, layout=name_button_layout)
        param_name17.style.button_color = 'lightgreen'

        self.number_of_scouts = IntText(
          value=10,
          step=1,
          style=style, layout=widget_layout)

        units_button1 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button1.style.button_color = 'lightgreen'
        units_button2 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        units_button2.style.button_color = 'tan'
        units_button3 = Button(description='min', disabled=True, layout=units_button_layout) 
        units_button3.style.button_color = 'lightgreen'
        units_button4 = Button(description='micron/min', disabled=True, layout=units_button_layout) 
        units_button4.style.button_color = 'tan'
        units_button5 = Button(description='micron/min', disabled=True, layout=units_button_layout) 
        units_button5.style.button_color = 'lightgreen'
        units_button6 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        units_button6.style.button_color = 'tan'
        units_button7 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button7.style.button_color = 'lightgreen'
        units_button8 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        units_button8.style.button_color = 'tan'
        units_button9 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        units_button9.style.button_color = 'lightgreen'
        units_button10 = Button(description='min', disabled=True, layout=units_button_layout) 
        units_button10.style.button_color = 'tan'
        units_button11 = Button(description='micron/min', disabled=True, layout=units_button_layout) 
        units_button11.style.button_color = 'lightgreen'
        units_button12 = Button(description='micron/min', disabled=True, layout=units_button_layout) 
        units_button12.style.button_color = 'tan'
        units_button13 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button13.style.button_color = 'lightgreen'
        units_button14 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button14.style.button_color = 'tan'
        units_button15 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button15.style.button_color = 'lightgreen'
        units_button16 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button16.style.button_color = 'tan'
        units_button17 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button17.style.button_color = 'lightgreen'

        desc_button1 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button1.style.button_color = 'lightgreen'
        desc_button2 = Button(description='the oxygen uptake rate of infection cell', disabled=True, layout=desc_button_layout) 
        desc_button2.style.button_color = 'tan'
        desc_button3 = Button(description='the persistence time of infection cell motility', disabled=True, layout=desc_button_layout) 
        desc_button3.style.button_color = 'lightgreen'
        desc_button4 = Button(description='the migration speed of infection cell motility', disabled=True, layout=desc_button_layout) 
        desc_button4.style.button_color = 'tan'
        desc_button5 = Button(description='the migration bias of unattached infection cell', disabled=True, layout=desc_button_layout) 
        desc_button5.style.button_color = 'lightgreen'
        desc_button6 = Button(description='the apoptosis rate of infection cell', disabled=True, layout=desc_button_layout) 
        desc_button6.style.button_color = 'tan'
        desc_button7 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button7.style.button_color = 'lightgreen'
        desc_button8 = Button(description='the oxygen uptake rate of scout cell', disabled=True, layout=desc_button_layout) 
        desc_button8.style.button_color = 'tan'
        desc_button9 = Button(description='the apoptosis rate of scout cell', disabled=True, layout=desc_button_layout) 
        desc_button9.style.button_color = 'lightgreen'
        desc_button10 = Button(description='the persistence time of scout cell motility', disabled=True, layout=desc_button_layout) 
        desc_button10.style.button_color = 'tan'
        desc_button11 = Button(description='the migration speed of scout cell motility', disabled=True, layout=desc_button_layout) 
        desc_button11.style.button_color = 'lightgreen'
        desc_button12 = Button(description='the migration bias of unattached scout cell', disabled=True, layout=desc_button_layout) 
        desc_button12.style.button_color = 'tan'
        desc_button13 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button13.style.button_color = 'lightgreen'
        desc_button14 = Button(description='the relative adhesion of scout cell', disabled=True, layout=desc_button_layout) 
        desc_button14.style.button_color = 'tan'
        desc_button15 = Button(description='the relative repulsion of scout cell', disabled=True, layout=desc_button_layout) 
        desc_button15.style.button_color = 'lightgreen'
        desc_button16 = Button(description='the number of infection cell', disabled=True, layout=desc_button_layout) 
        desc_button16.style.button_color = 'tan'
        desc_button17 = Button(description='the number of scouts cell', disabled=True, layout=desc_button_layout) 
        desc_button17.style.button_color = 'lightgreen'

        row1 = [param_name1, self.random_seed, units_button1, desc_button1] 
        row2 = [param_name2, self.infection_o2_relative_uptake, units_button2, desc_button2] 
        row3 = [param_name3, self.infection_persistence_time, units_button3, desc_button3] 
        row4 = [param_name4, self.infection_migration_speed, units_button4, desc_button4] 
        row5 = [param_name5, self.unattached_infection_migration_bias, units_button5, desc_button5] 
        row6 = [param_name6, self.infection_apoptosis_rate, units_button6, desc_button6] 
        row7 = [param_name7, self.infection_relative_cycle_entry_rate, units_button7, desc_button7] 
        row8 = [param_name8, self.scout_o2_relative_uptake, units_button8, desc_button8] 
        row9 = [param_name9, self.scout_apoptosis_rate, units_button9, desc_button9] 
        row10 = [param_name10, self.scout_persistence_time, units_button10, desc_button10] 
        row11 = [param_name11, self.scout_migration_speed, units_button11, desc_button11] 
        row12 = [param_name12, self.unattached_scout_migration_bias, units_button12, desc_button12] 
        row13 = [param_name13, self.scout_relative_cycle_entry_rate, units_button13, desc_button13] 
        row14 = [param_name14, self.scout_relative_adhesion, units_button14, desc_button14] 
        row15 = [param_name15, self.scout_relative_repulsion, units_button15, desc_button15] 
        row16 = [param_name16, self.number_of_infection, units_button16, desc_button16] 
        row17 = [param_name17, self.number_of_scouts, units_button17, desc_button17] 

        box_layout = Layout(display='flex', flex_flow='row', align_items='stretch', width='100%')
        box1 = Box(children=row1, layout=box_layout)
        box2 = Box(children=row2, layout=box_layout)
        box3 = Box(children=row3, layout=box_layout)
        box4 = Box(children=row4, layout=box_layout)
        box5 = Box(children=row5, layout=box_layout)
        box6 = Box(children=row6, layout=box_layout)
        box7 = Box(children=row7, layout=box_layout)
        box8 = Box(children=row8, layout=box_layout)
        box9 = Box(children=row9, layout=box_layout)
        box10 = Box(children=row10, layout=box_layout)
        box11 = Box(children=row11, layout=box_layout)
        box12 = Box(children=row12, layout=box_layout)
        box13 = Box(children=row13, layout=box_layout)
        box14 = Box(children=row14, layout=box_layout)
        box15 = Box(children=row15, layout=box_layout)
        box16 = Box(children=row16, layout=box_layout)
        box17 = Box(children=row17, layout=box_layout)

        self.tab = VBox([
          box1,
          box2,
          box3,
          box4,
          box5,
          box6,
          box7,
          box8,
          box9,
          box10,
          box11,
          box12,
          box13,
          box14,
          box15,
          box16,
          box17,
        ])

    # Populate the GUI widgets with values from the XML
    def fill_gui(self, xml_root):
        uep = xml_root.find('.//user_parameters')  # find unique entry point into XML
        self.random_seed.value = int(uep.find('.//random_seed').text)
        self.infection_o2_relative_uptake.value = float(uep.find('.//infection_o2_relative_uptake').text)
        self.infection_persistence_time.value = float(uep.find('.//infection_persistence_time').text)
        self.infection_migration_speed.value = float(uep.find('.//infection_migration_speed').text)
        self.unattached_infection_migration_bias.value = float(uep.find('.//unattached_infection_migration_bias').text)
        self.infection_apoptosis_rate.value = float(uep.find('.//infection_apoptosis_rate').text)
        self.infection_relative_cycle_entry_rate.value = float(uep.find('.//infection_relative_cycle_entry_rate').text)
        self.scout_o2_relative_uptake.value = float(uep.find('.//scout_o2_relative_uptake').text)
        self.scout_apoptosis_rate.value = float(uep.find('.//scout_apoptosis_rate').text)
        self.scout_persistence_time.value = float(uep.find('.//scout_persistence_time').text)
        self.scout_migration_speed.value = float(uep.find('.//scout_migration_speed').text)
        self.unattached_scout_migration_bias.value = float(uep.find('.//unattached_scout_migration_bias').text)
        self.scout_relative_cycle_entry_rate.value = float(uep.find('.//scout_relative_cycle_entry_rate').text)
        self.scout_relative_adhesion.value = float(uep.find('.//scout_relative_adhesion').text)
        self.scout_relative_repulsion.value = float(uep.find('.//scout_relative_repulsion').text)
        self.number_of_infection.value = int(uep.find('.//number_of_infection').text)
        self.number_of_scouts.value = int(uep.find('.//number_of_scouts').text)


    # Read values from the GUI widgets to enable editing XML
    def fill_xml(self, xml_root):
        uep = xml_root.find('.//user_parameters')  # find unique entry point into XML 
        uep.find('.//random_seed').text = str(self.random_seed.value)
        uep.find('.//infection_o2_relative_uptake').text = str(self.infection_o2_relative_uptake.value)
        uep.find('.//infection_persistence_time').text = str(self.infection_persistence_time.value)
        uep.find('.//infection_migration_speed').text = str(self.infection_migration_speed.value)
        uep.find('.//unattached_infection_migration_bias').text = str(self.unattached_infection_migration_bias.value)
        uep.find('.//infection_apoptosis_rate').text = str(self.infection_apoptosis_rate.value)
        uep.find('.//infection_relative_cycle_entry_rate').text = str(self.infection_relative_cycle_entry_rate.value)
        uep.find('.//scout_o2_relative_uptake').text = str(self.scout_o2_relative_uptake.value)
        uep.find('.//scout_apoptosis_rate').text = str(self.scout_apoptosis_rate.value)
        uep.find('.//scout_persistence_time').text = str(self.scout_persistence_time.value)
        uep.find('.//scout_migration_speed').text = str(self.scout_migration_speed.value)
        uep.find('.//unattached_scout_migration_bias').text = str(self.unattached_scout_migration_bias.value)
        uep.find('.//scout_relative_cycle_entry_rate').text = str(self.scout_relative_cycle_entry_rate.value)
        uep.find('.//scout_relative_adhesion').text = str(self.scout_relative_adhesion.value)
        uep.find('.//scout_relative_repulsion').text = str(self.scout_relative_repulsion.value)
        uep.find('.//number_of_infection').text = str(self.number_of_infection.value)
        uep.find('.//number_of_scouts').text = str(self.number_of_scouts.value)
