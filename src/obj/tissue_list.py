import csv
# Define the CustomObject class
class CustomObject:
    def __init__(self, name):
        self.name = name
        self.values = []

    def add_value(self, value):
        self.values.append(value)

class TissueList:
    def __init__(self):
        # Read the CSV file
        with open('src/data_input/input_CCLE.csv', 'r') as csvfile:
            reader = csv.reader(csvfile)
            csv_data = list(reader)

        # Dictionary to store objects
        self.objects_dict = {}
        self.elements_list = []

        # Process each row
        for row in csv_data:
            # Skip empty rows
            if not row:
                continue

            # Split the first element of the row at the first underscore
            name, elements = row[0].split('_', 1)
            self.elements_list.append(elements)

        tissue_list = list(set(self.elements_list))

        for name in tissue_list:
            custom_object = CustomObject(name)
            self.objects_dict[name] = custom_object

        for row in csv_data:
            # Skip empty rows
            if not row:
                continue

            # Split the first element of the row at the first underscore
            name, elements = row[0].split('_', 1)
            self.objects_dict[elements].add_value(name)

    def get_tissue_list(self):
        return sorted([obj for obj in self.objects_dict])

    def get_ccle_list(self, tissue_name):
        return sorted(self.objects_dict[tissue_name].values)

tissue_list = TissueList()
#print ('stampo la lista dei tessuti:')
#print(tissue_list.get_tissue_list())
#print('-----------------')
#print ('stampo la lista delle linee cellulari per CENTRAL_NERVOUS_SYSTEM:')
#print(tissue_list.get_ccle_list('CENTRAL_NERVOUS_SYSTEM'))
