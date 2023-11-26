import os

import rpy2.robjects as robjects
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QFont
from PyQt5.QtWidgets import (QWidget, QPushButton, QVBoxLayout,
                             QLabel, QComboBox, QRadioButton, QGroupBox, QHBoxLayout, QTextEdit,
                             QFrame, QMessageBox, QToolTip)
from rpy2.robjects import pandas2ri
from tabulate import tabulate

from src.obj.tissue_list import TissueList
from src.window.dataframe_window import DataFrameViewer


class RepurposingWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.element_combo = None
        self.tissue_combo = None
        self.selected_method_text = None
        self.timeline_text = None
        self.tissue_list = TissueList()
        self.viewer = None
        self.init_ui()

    def update_secondary_combo(self):
        # Update the first combo box based on the selected tissue
        selected_tissue = self.tissue_combo.currentText()
        self.element_combo.clear()
        self.element_combo.addItems(self.tissue_list.get_ccle_list(selected_tissue))

    def handle_selection_change(self, index, console_output):
        console_output.clear()
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()

        # You can also access the console_output if needed
        console_output.append(f"--------------")
        console_output.append(f"Selected Tissue: {selected_tissue}")
        console_output.append(f"Selected CCLE: {selected_element}")
        console_output.append(f"--------------")

    def handle_mouse_hover(self, event):
        # Show the tooltip when the mouse hovers over the label
        QToolTip.showText(self.mapToGlobal(event.pos()), 'MSSG (Most Statistically Significative Genes)')

    def init_ui(self):

        # Add a logo at the top
        logo_label = QLabel(self)
        pixmap = QPixmap('image/logo.png')  # Replace 'path_to_your_logo.png' with the actual path to your logo image
        scaled_pixmap = pixmap.scaledToWidth(100)  # Adjust the width as needed
        logo_label.setPixmap(scaled_pixmap)
        logo_label.setAlignment(Qt.AlignLeft)

        intestazione_label = QLabel("Drug Repurposing from existing Experiment")
        ifont = QFont("Arial", 14)
        ifont.setBold(True)
        intestazione_label.setFont(ifont)
        intestazione_label.setAlignment(Qt.AlignCenter | Qt.AlignTop)

        # Add a horizontal separator
        separator_1 = QFrame()
        separator_1.setFrameShape(QFrame.HLine)
        separator_1.setFrameShadow(QFrame.Sunken)

        separator_2 = QFrame()
        separator_2.setFrameShape(QFrame.HLine)
        separator_2.setFrameShadow(QFrame.Sunken)

        separator_3 = QFrame()
        separator_3.setFrameShape(QFrame.HLine)
        separator_3.setFrameShadow(QFrame.Sunken)

        separator_4 = QFrame()
        separator_4.setFrameShape(QFrame.HLine)
        separator_4.setFrameShadow(QFrame.Sunken)

        separator_4 = QFrame()
        separator_4.setFrameShape(QFrame.HLine)
        separator_4.setFrameShadow(QFrame.Sunken)

        separator_5 = QFrame()
        separator_5.setFrameShape(QFrame.HLine)
        separator_5.setFrameShadow(QFrame.Sunken)

        separator_6 = QFrame()
        separator_6.setFrameShape(QFrame.HLine)
        separator_6.setFrameShadow(QFrame.Sunken)

        # Create a combo box for TESSUTO
        tissue_label = QLabel('Choose a Tissue:')
        font = QFont("Arial", 12)
        font.setBold(True)
        tissue_label.setFont(font)
        self.tissue_combo = QComboBox(self)
        # Ottengo la lista dei Tessuti
        self.tissue_combo.addItems(self.tissue_list.get_tissue_list())
        predefined_tissue = "LUNG"
        index = self.tissue_combo.findText(predefined_tissue)
        if index != -1:
            # If the item is found, set it as the current index
            self.tissue_combo.setCurrentIndex(index)
        #self.tissue_combo.setDisabled(True)
        self.tissue_combo.currentIndexChanged.connect(self.update_secondary_combo)

        # Create a QTextEdit for console output
        self.console_output = QTextEdit(self)
        self.console_output.setReadOnly(True)  # Set it to read-only

        # Create the second combo box for selecting the element with the chosen tissue
        element_label = QLabel("Choose a Cellular Line...")
        element_label.setFont(font)
        self.element_combo = QComboBox(self)
        predefined_element = "A549"
        self.element_combo.addItem(predefined_element)
        self.element_combo.setCurrentText(predefined_element)
        #self.element_combo.setDisabled(True)
        self.console_output.append(f"--------------")
        self.console_output.append(f"Selected Tissue: {predefined_tissue}")
        self.console_output.append(f"Selected CCLE: {predefined_element}")
        self.console_output.append(f"--------------")

        # Connect the clicked signal of the QListView to a custom function
        self.tissue_combo.currentIndexChanged.connect(
            lambda index: self.handle_selection_change(index, self.console_output))
        self.element_combo.currentIndexChanged.connect(
            lambda index: self.handle_selection_change(index, self.console_output))

        # Create radio buttons for METODO
        metodo_label = QLabel('Choose a method ')
        metodo_label.setToolTip('MSSG (Most Statistically Significative Genes)')
        metodo_label.mousePressEvent = self.handle_mouse_hover

        metodo_label.setFont(font)
        self.metodo_group = QGroupBox()
        metodo_layout = QHBoxLayout()
        self.metodo_50mssg = QRadioButton('50 Most MSSG')
        self.metodo_100mssg = QRadioButton('100 Most MSSG')
        self.metodo_150mssg = QRadioButton('150 Most MSSG')
        self.metodo_BCmssg = QRadioButton('Bin Chen Method')
        metodo_layout.addWidget(self.metodo_50mssg)
        metodo_layout.addWidget(self.metodo_100mssg)
        metodo_layout.addWidget(self.metodo_150mssg)
        metodo_layout.addWidget(self.metodo_BCmssg)
        self.metodo_group.setLayout(metodo_layout)

        # Create radio buttons for HOURS
        time_esperiment_label = QLabel('Experiment duration:')
        time_esperiment_label.setFont(font)
        self.time_esperiment_group = QGroupBox()
        time_esperiment_layout = QHBoxLayout()
        self.time_esperiment_6h = QRadioButton('6h')
        self.time_esperiment_24h = QRadioButton('24h')
        self.time_esperiment_6h_24h = QRadioButton('6h_24h')
        self.time_esperiment_6h_24h_c = QRadioButton('6h_24_corrected')
        time_esperiment_layout.addWidget(self.time_esperiment_6h)
        time_esperiment_layout.addWidget(self.time_esperiment_24h)
        time_esperiment_layout.addWidget(self.time_esperiment_6h_24h)
        time_esperiment_layout.addWidget(self.time_esperiment_6h_24h_c)
        self.time_esperiment_group.setLayout(time_esperiment_layout)

        # Create a button
        self.button = QPushButton('Extract data', self)

        # Create a vertical layout for the lateral panels
        # lateral_layout = QVBoxLayout()

        # Create a vertical layout for the main content
        main_layout = QVBoxLayout()
        logo_layout = QHBoxLayout()
        logo_layout.addWidget(logo_label)
        logo_layout.addWidget(intestazione_label)
        main_layout.addLayout(logo_layout)
        # inserisco le combobox per la selezione del tessuto e della linea cellulare
        main_layout.addWidget(separator_1)
        main_layout.addWidget(tissue_label)
        main_layout.addWidget(self.tissue_combo)
        main_layout.addWidget(element_label)
        main_layout.addWidget(self.element_combo)

        main_layout.addWidget(separator_3)
        main_layout.addWidget(metodo_label)
        main_layout.addWidget(self.metodo_group)
        main_layout.addWidget(separator_4)
        main_layout.addWidget(time_esperiment_label)
        main_layout.addWidget(self.time_esperiment_group)
        main_layout.addWidget(separator_5)
        main_layout.addWidget(self.button)
        main_layout.addWidget(separator_6)
        main_layout.addWidget(self.console_output)

        # Create a horizontal layout to combine the lateral and main layouts
        combined_layout = QHBoxLayout()
        # combined_layout.addLayout(lateral_layout)
        combined_layout.addLayout(main_layout)

        # Set the combined layout for the main window
        self.setLayout(combined_layout)

        # Connect the button click event to a function
        self.button.clicked.connect(lambda: self.on_button_click(self.console_output))

        self.metodo_50mssg.clicked.connect(
            lambda checked: self.on_method_selected_in_repurposing('50 Most MSSG', self.console_output))
        self.metodo_100mssg.clicked.connect(
            lambda checked: self.on_method_selected_in_repurposing('100 Most MSSG', self.console_output))
        self.metodo_150mssg.clicked.connect(
            lambda checked: self.on_method_selected_in_repurposing('150 Most MSSG', self.console_output))
        self.metodo_BCmssg.clicked.connect(
            lambda checked: self.on_method_selected_in_repurposing('BinChen selected', self.console_output))

        self.time_esperiment_6h.clicked.connect(
            lambda checked: self.on_timeline_selected_in_repurposing('6h', self.console_output))
        self.time_esperiment_24h.clicked.connect(
            lambda checked: self.on_timeline_selected_in_repurposing('24h', self.console_output))
        self.time_esperiment_6h_24h.clicked.connect(
            lambda checked: self.on_timeline_selected_in_repurposing('6h_24h', self.console_output))
        self.time_esperiment_6h_24h_c.clicked.connect(
            lambda checked: self.on_timeline_selected_in_repurposing('6h_24h_c', self.console_output))

    def on_method_selected_in_repurposing(self, checked, console_output):
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()
        self.selected_method_text = checked
        console_output.clear()
        console_output.append(f"--------------")
        console_output.append(f"Selected Tissue: {selected_tissue}")
        console_output.append(f"Selected CCLE: {selected_element}")
        console_output.append(f"Selected Method: {self.selected_method_text}")
        console_output.append(f"--------------")

    def on_timeline_selected_in_repurposing(self, checked_time, console_output):
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()
        self.timeline_text = checked_time
        console_output.clear()
        console_output.append(f"--------------")
        console_output.append(f"Selected Tissue: {selected_tissue}")
        console_output.append(f"Selected CCLE: {selected_element}")
        console_output.append(f"Selected Method: {self.selected_method_text}")
        console_output.append(f"Experiment duration: {self.timeline_text}")
        console_output.append(f"--------------")

    def on_button_click(self, console_output):
        m = self.selected_method_text
        t = self.timeline_text
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()
        scenario = ''
        console_output.append("Processing...")
        if m == '50 Most MSSG':
            if t == '6h':
                console_output.append(f"Scenario 1 : {m},{t}")
                scenario = 'scenario1'
            elif t == '24h':
                console_output.append(f"Scenario 2 : {m},{t}")
                scenario = 'scenario2'
            elif t == '6h_24h':
                console_output.append(f"Scenario 3 : {m},{t}")
                scenario = 'scenario3'
            elif t == '6h_24h_c':
                console_output.append(f"Scenario 4 : {m},{t}")
                scenario = 'scenario4'
        elif m == '100 Most MSSG':
            if t == '6h':
                console_output.append(f"Scenario 5 : {m},{t}")
                scenario = 'scenario5'
            elif t == '24h':
                console_output.append(f"Scenario 6 : {m},{t}")
                scenario = 'scenario6'
            elif t == '6h_24h':
                console_output.append(f"Scenario 7 : {m},{t}")
                scenario = 'scenario7'
            elif t == '6h_24h_c':
                console_output.append(f"Scenario 8 : {m},{t}")
                scenario = 'scenario8'
        elif m == '150 Most MSSG':
            if t == '6h':
                console_output.append(f"Scenario 9 : {m},{t}")
                scenario = 'scenario9'
            elif t == '24h':
                console_output.append(f"Scenario 10 : {m},{t}")
                scenario = 'scenario10'
            elif t == '6h_24h':
                console_output.append(f"Scenario 11 : {m},{t}")
                scenario = 'scenario11'
            elif t == '6h_24h_c':
                console_output.append(f"Scenario 12 : {m},{t}")
                scenario = 'scenario12'
        elif m == 'BinChen selected':
            if t == '6h':
                console_output.append(f"Scenario 13 : {m},{t}")
                scenario = 'scenario13'
            elif t == '24h':
                console_output.append(f"Scenario 14 : {m},{t}")
                scenario = 'scenario14'
            elif t == '6h_24h':
                console_output.append(f"Scenario 15 : {m},{t}")
                scenario = 'scenario15'
            elif t == '6h_24h_c':
                console_output.append(f"Scenario 16 : {m},{t}")
                scenario = 'scenario16'

        # Load the RDS file
        rds_file_path = ('src/r_code/Output/4.AUC_versus_connectivity_score/'+
                         selected_element +'_'+
                         selected_tissue+'/'+
                         scenario + '.rds')
        if not os.path.exists(rds_file_path):
            console_output.append(f" Path not found: {rds_file_path}")
            console_output.append(f" Seems like there's no experiment for {selected_element}_{selected_tissue}")
            return
        readRDS = robjects.r['readRDS']
        df = readRDS(rds_file_path)
        df = pandas2ri.rpy2py_dataframe(df)

        # do something with the dataframe
        # Filter rows based on the condition
        # filtered_df = df[df['connectivity_score_with_LINCS_profile.pvalue'].between(0.0001, 0.005)]

        # Sort the filtered DataFrame
        # FIltro per valori P.Value p<0,05 & p>0,0001
        # sorted_df = filtered_df.sort_values(by=['connectivity_score_with_LINCS_profile.pvalue'])
        df = df[['drug_name',
                'connectivity_score_with_LINCS_profile',
                'connectivity_score_with_LINCS_profile.pvalue',
                'alt_connectivity_score_with_LINCS_profile',
                'alt_connectivity_score_with_LINCS_profile.pvalue']]

        df = df.sort_values(['connectivity_score_with_LINCS_profile','connectivity_score_with_LINCS_profile.pvalue'], ascending=[True, True])

        df_top10 = df[['drug_name', 'connectivity_score_with_LINCS_profile']].head(10)
        # print(df[['drug_name', 'connectivity_score_with_LINCS_profile']].head(10))

        print(tabulate(df_top10, headers='keys', tablefmt='psql'))
        console_output.append(f" Selected drugs \n: {tabulate(df_top10, headers='keys', tablefmt='psql')}")

        selected_data = self.get_selected_data()
        self.viewer = DataFrameViewer(df, selected_data)
        self.viewer.show()

        # self.show_message_box(scenario)
        # self.open_second_window()


    def get_selected_data(self):
        # Retrieve data from your widgets and return it
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()
        selected_method = self.selected_method_text  # Assuming you've set this somewhere
        selected_timeline = self.timeline_text  # Assuming you've set this somewhere

        # Combine the data into a dictionary or whatever format you prefer
        selected_data = {
            "selected_tissue": selected_tissue,
            "selected_element": selected_element,
            "selected_method": selected_method,
            "selected_timeline": selected_timeline
            # Add more keys as needed
        }

        return selected_data
