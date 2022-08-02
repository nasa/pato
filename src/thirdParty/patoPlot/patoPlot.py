from pylab import *

# pato plot class reads the list of data files
class patoPlot:
    def __init__(self, fileNameList):
        self.fileNameList = fileNameList # list of data file
        self.dataList = [] # list of data
        self.titleList = [] # list of titles
        # reads the list of data files
        for fileName in self.fileNameList: 
            data,title= self.readDataPATO(fileName) # reads data file
            self.dataList.append(data) # store data
            self.titleList.append(title) # store title

    # return the number in a string
    def num(self, s):
        table=[] # creates table
        try:
            for i in s:
                table.append(float(i)) # add the number i from string s in the table
            return table
        except ValueError:
            return table # error if there is no number in the string s

    # read data from filename
    def readDataPATO(self, fileName):
        data = [] # data
        title = [] # title
        with open(fileName) as f:
            content = f.readlines() # read lines
            i = 0
            for line in content:
                if i == 0:
                    title=line.split() # take title
                if i > 0:
                    values = line.split()
                    data.append(self.num(values)) # take values from string
                i = i + 1
        data=list(filter(None,data))
        first = len(data[0]) # length of the first column
        addColumn = 1
        for row in data:
            lenRow = len(row) # length of the row
            if lenRow<first:
                for i in range(0,first-lenRow):
                    row.insert(i+addColumn,None) # add None if the length is different 

        data=array(data)
        return data,title

