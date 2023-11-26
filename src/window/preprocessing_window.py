import subprocess
import traceback

import yaml
from PyQt5.QtCore import *
from PyQt5.QtCore import Qt, QRunnable, QObject
from PyQt5.QtGui import QPixmap, QFont
from PyQt5.QtWidgets import (QWidget, QPushButton, QVBoxLayout,
                             QLabel, QComboBox, QLineEdit, QRadioButton, QGroupBox, QHBoxLayout, QTextEdit,
                             QFrame, QToolTip)

from src.obj.tissue_list import TissueList


class PreProcessingWindow(QWidget):
    def __init__(self, *args, **kwargs):
        super(PreProcessingWindow, self).__init__(*args, **kwargs)
        # ... (rest of the initialization)

        self.threadpool = QThreadPool()
        self.threadpool.setMaxThreadCount(1)
        print("Multithreading with maximum %d threads" % self.threadpool.maxThreadCount())

        self.selected_sample = None
        self.sample_input = None
        self.element_combo = None
        self.tissue_combo = None
        self.selected_method_text = None
        self.timeline_text = None
        self.tissue_list = TissueList()
        self.init_ui()

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

        intestazione_label = QLabel("Preprocessing Tissue with specific CCLE")
        ifont = QFont("Arial", 14)
        ifont.setBold(True)
        intestazione_label.setFont(ifont)
        intestazione_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

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
        self.tissue_combo.currentIndexChanged.connect(self.update_secondary_combo)

        # Create the second combo box for selecting the element with the chosen tissue
        element_label = QLabel('Choose a Cellular Line:')
        element_label.setFont(font)
        self.element_combo = QComboBox(self)
        self.element_combo.addItem("Choose a Cellular Line...")


        # Create an input text field for Sample
        sample_label = QLabel('Insert samples (preprocessing):')

        sample_label.setFont(font)
        self.sample_input = QLineEdit(self)

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
        self.metodo_group.setDisabled(True)
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
        self.time_esperiment_group.setDisabled(True)

        self.time_esperiment_group.setLayout(time_esperiment_layout)

        # Create a button
        self.button = QPushButton('Preprocess', self)

        # Create a QTextEdit for console output
        self.console_output = QTextEdit(self)
        self.console_output.setReadOnly(True)  # Set it to read-only

        # Create a vertical layout for the lateral panels
        # lateral_layout = QVBoxLayout()

        # Create a vertical layout for the main content
        main_layout = QVBoxLayout()
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
        main_layout.addWidget(separator_2)
        main_layout.addWidget(sample_label)
        main_layout.addWidget(self.sample_input)
        main_layout.addWidget(separator_3)
        #Disabilito perchè non serve

        #main_layout.addWidget(metodo_label)
        #main_layout.addWidget(self.metodo_group)
        #main_layout.addWidget(separator_4)
        #Disabilito perchè non serve
        #main_layout.addWidget(time_esperiment_label)
        #main_layout.addWidget(self.time_esperiment_group)
        #main_layout.addWidget(separator_5)
        main_layout.addWidget(self.button)
        main_layout.addWidget(separator_6)
        main_layout.addWidget(self.console_output)

        # Create a horizontal layout to combine the lateral and main layouts
        combined_layout = QHBoxLayout()
        # combined_layout.addLayout(lateral_layout)
        combined_layout.addLayout(main_layout)

        # Set the combined layout for the main window
        self.setLayout(combined_layout)

        # Connect the clicked signal of the QListView to a custom function
        self.tissue_combo.currentIndexChanged.connect(
            lambda index: self.on_tissue_selected_in_preprop(index, self.console_output))

        self.element_combo.currentIndexChanged.connect(
            lambda index: self.on_tissue_selected_in_preprop(index, self.console_output))

        # Connect the button click event to a function
        self.button.clicked.connect(lambda: self.on_button_click_in_preprop(self.console_output))


        self.sample_input.textChanged.connect(
            lambda: self.on_sample_selected_in_preprop(self.console_output))
        self.sample_input.returnPressed.connect(
            lambda: self.on_sample_selected_in_preprop(self.console_output))


        self.metodo_50mssg.clicked.connect(
            lambda checked: self.on_method_selected_in_preprop('50 Most MSSG', self.console_output))
        self.metodo_100mssg.clicked.connect(
            lambda checked: self.on_method_selected_in_preprop('100 Most MSSG', self.console_output))
        self.metodo_150mssg.clicked.connect(
            lambda checked: self.on_method_selected_in_preprop('150 Most MSSG', self.console_output))
        self.metodo_BCmssg.clicked.connect(
            lambda checked: self.on_method_selected_in_preprop('BinChen selected', self.console_output))

        self.time_esperiment_6h.clicked.connect(
            lambda checked: self.on_timeline_selected_in_preprop('6h', self.console_output))
        self.time_esperiment_24h.clicked.connect(
            lambda checked: self.on_timeline_selected_in_preprop('24h', self.console_output))
        self.time_esperiment_6h_24h.clicked.connect(
            lambda checked: self.on_timeline_selected_in_preprop('6h_24h', self.console_output))
        self.time_esperiment_6h_24h_c.clicked.connect(
            lambda checked: self.on_timeline_selected_in_preprop('6h_24h_c', self.console_output))

    def update_secondary_combo(self):
        # Update the first combo box based on the selected tissue
        selected_tissue = self.tissue_combo.currentText()
        self.element_combo.clear()
        self.element_combo.addItems(self.tissue_list.get_ccle_list(selected_tissue))

    def on_tissue_selected_in_preprop(self, index, console_output):
        console_output.clear()
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()

        # You can also access the console_output if needed
        console_output.append(f"--------------")
        console_output.append(f"Selected Tissue: {selected_tissue}")
        console_output.append(f"Selected CCLE: {selected_element}")
        console_output.append(f"--------------")
    def on_sample_selected_in_preprop(self, console_output):
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()
        self.selected_sample = self.sample_input.text()
        console_output.clear()
        console_output.append(f"--------------")
        console_output.append(f"Selected Tissue: {selected_tissue}")
        console_output.append(f"Selected CCLE: {selected_element}")
        console_output.append(f"Samples: {self.selected_sample}")
        console_output.append(f"--------------")

    def on_method_selected_in_preprop(self, checked, console_output):
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()
        self.selected_method_text = checked
        console_output.clear()
        console_output.append(f"--------------")
        console_output.append(f"Selected Tissue: {selected_tissue}")
        console_output.append(f"Selected CCLE: {selected_element}")
        console_output.append(f"Samples: {self.selected_sample}")
        console_output.append(f"Selected Method: {self.selected_method_text}")
        console_output.append(f"--------------")

    def on_timeline_selected_in_preprop(self, checked_time, console_output):
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()
        self.timeline_text = checked_time
        console_output.clear()
        console_output.append(f"--------------")
        console_output.append(f"Selected Tissue: {selected_tissue}")
        console_output.append(f"Selected CCLE: {selected_element}")
        console_output.append(f"Samples: {self.selected_sample}")
        console_output.append(f"Selected Method: {self.selected_method_text}")
        console_output.append(f"Experiment Duration: {self.timeline_text}")
        console_output.append(f"--------------")

    def on_button_click_in_preprop(self, console_output):
        selected_tissue = self.tissue_combo.currentText()
        selected_element = self.element_combo.currentText()

        if selected_element == "Choose a Cellular Line...":
            console_output.clear()
            console_output.append('Choose a Cellular Line!')
            self.tissue_combo.setFocus()
            return

        if self.selected_sample == None:
            console_output.clear()
            console_output.append('Select samples!')
            return

        console_output.clear()
        console_output.append('Selected Configuration:')
        console_output.append(f"--------------")
        console_output.append(f"Selected Tissue: {selected_tissue}")
        console_output.append(f"Selected CCLE: {selected_element}")
        console_output.append(f"Samples: {self.selected_sample}")
        console_output.append(f"--------------")

        path = 'src/r_code/test_python/'
        r_script = 'main.R'


        script_arguments = [str(selected_tissue), #args[1]
                            str(selected_element),#args[2]
                            str(self.selected_sample),#args[3]
                            str(path)#args[4]
                            ]
        #Creazioen file YAML
        linea_ccle = selected_element+'_'+selected_tissue

        parametri = {
            'linea_ccle' :
                linea_ccle,
            'campioni' :
                self.selected_sample
        }
        percorso_file_config = 'src/config.yaml'
        # Scrivi i parametri nel file config.yaml utilizzando PyYAML
        with open(percorso_file_config, 'w') as file_config:
            yaml.dump(parametri, file_config, default_flow_style=False)

        print(f"File config.yaml created in {percorso_file_config}")

        # Specifica il percorso del tuo Snakefile e del file config.yaml
        percorso_snakefile = 'src/snakefile'

        # Costruisci il comando per eseguire Snakemake
        comando_snakemake = f"snakemake --cores 1 --snakefile  {percorso_snakefile} --configfile {percorso_file_config}"

        # Esegui il comando utilizzando subprocess
        try:
            subprocess.run(comando_snakemake, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error during execution of Snakemake: {e}")


        #print(r_script)
        #worker = Worker(self.run_script_R, path + r_script, script_arguments, console_output)
        #worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
        #worker.signals.progress.connect(self.progress_fn)
        #self.threadpool.start(worker)

    def thread_complete(self):
        print("THREAD COMPLETE!")

    def progress_fn(self):
        print("Connessione...")

    def print_output(self):
        print('Thread Concluso!')

    def run_script_R(self, r_script_path, r_script_argument,console_output):
        try:
            # Run the R script using the R interpreter
            result = subprocess.run(['Rscript', r_script_path] + r_script_argument, check=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True)
            output = ''
            if result.returncode == 0:
                output = result.stdout
                print(f"Output:\n{output}")
            else:
                error = result.stderr
                print(f"Error:\n{error}")
            console_output.append(r_script_path+"R script executed successfully.\n" + output)

        except subprocess.CalledProcessError as e:
            console_output.append(f"Error executing R script:\n{e.stderr}")
        except Exception as e:
            console_output.append(f"An error occurred: {str(e)}")

class Worker(QRunnable):
    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
            self.signals.result.emit(result)
        except Exception as e:
            self.signals.error.emit((type(e), str(e), traceback.format_exc()))
        finally:
            self.signals.finished.emit()

class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    progress
        int indicating % progress

    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)

