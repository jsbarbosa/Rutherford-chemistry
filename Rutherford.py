import numpy as numpy
import matplotlib.pyplot as plt
import os
import wx
import sys
#from mpl_toolkits.mplot3d import *
from threading import Thread
#from wx.lib.pubsub import pub
import time

from scipy.signal import find_peaks_cwt
from scipy.signal import wiener

wildcard = "Semicolon-separated values (*.csv)|*.csv|""Text (*.txt)|*.txt|" \
            "All files (*.*)|*.*"
wildcard2 = "Text (*.txt)|*.txt|""Comma-separated values (*.csv)|*.csv|" \
            "All files (*.*)|*.*"
               
class initialRunThread(Thread):
    global mainWindow
    def __init__(self):
        Thread.__init__(self)
        self.daemon = True
        self.start()
    def run(self):
        mainWindow.initialRun()
        
class cargarThread(Thread):
    global mainWindow
    def __init__(self):
        Thread.__init__(self)
        self.daemon = True
        self.start()
    def run(self):
        mainWindow.cargar()
        
class massEspectrumForPeaksThread(Thread):
    global mainWindow
    def __init__(self):
        Thread.__init__(self)
        self.daemon = True
        self.start()
    def run(self):
        mainWindow.massEspectrumForPeaks()
        
class completeAnalysisThread(Thread):
    global mainWindow
    def __init__(self):
        Thread.__init__(self)
        self.start()
    def run(self):
        mainWindow.completeAnalysis()
        if mainWindow.darAutoClose():
            os._exit(1)
        print "Complete Analysis is over."
        
class windowClass(wx.Frame):
    def __init__(self, *args, **kwargs):
        no_resize = wx.DEFAULT_FRAME_STYLE & ~ (wx.RESIZE_BORDER | wx.RESIZE_BOX | wx.MAXIMIZE_BOX)
#        no_resize =  wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.STAY_ON_TOP
        super(windowClass, self).__init__(size = (700, 500), style= no_resize,  *args, **kwargs)
        self.userName = ""
        self.basicGUI()
        self.currentDirectory = os.getcwd()

    def basicGUI(self):
        nameBox = wx.TextEntryDialog(None, 'What is your name?', 'Welcome', 'UniAndes')
        if nameBox.ShowModal() == wx.ID_OK:
            self.userName = nameBox.GetValue()

        self.filePaths = []
        self.spectralPaths = []
        self.cromatogramas = []
        self.extention = ""
        
        ico = wx.Icon('icon.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(ico)

        panel = wx.Panel(self, size=(400,500))

        menuBar = wx.MenuBar()

        self.CreateStatusBar()

        #File initialization
        wx.StaticText(panel, label="Files(s)", pos=(230, 10))
        self.FileNameTextCtrl = wx.TextCtrl(panel, pos=(280, 8), size=(280, 20))
        self.openButton = wx.Button(panel, label="Browse", pos=(570, 5))
        wx.StaticText(panel, label="Output", pos=(230, 40))
        self.savePath = wx.TextCtrl(panel, value="", pos = (280,38), size = (280, 20))
        self.saveButton = wx.Button(panel, label="Browse", pos=(570, 35))
        self.automaticSaveCheckBox = wx.CheckBox(panel, label="Automatic saving", pos=(280, 60))

        #Compound Detection
        wx.StaticText(panel, label="Spectral Data", pos=(250, 230))
        self.spectralDataTextCtrl = wx.TextCtrl(panel, pos=(330, 227), size=(230, 20))
        self.spectralDataButton = wx.Button(panel, label="Browse", pos=(570, 225))
        self.spectralDataButtonRun = wx.Button(panel, label="Run", pos=(570, 255))

        #StaticBoxes
        simpleAnalysis = wx.StaticBox(panel, label="Simple Analysis", pos=(230, 90), size=(440, 110))
        compoundDetection = wx.StaticBox(panel, label="Compound Detection", pos=(230, 200), size=(440, 100))

        #Buttons
        cromatogramaButton = wx.Button(panel, label="Generate Chromatogram", pos=(250, 160), size = (150, 30))
        peaksMassSpectra = wx.Button(panel, label="Peaks Mass Spectra", pos=(470, 115), size = (150, 30))
        mass3dButton = wx.Button(panel, label="3D Mass Spectrum", pos=(470, 155), size = (150, 30))
        initialButtonRun = wx.Button(panel, label="Run", pos=(570, 65))
        simpleAnalysisButtonRun = wx.Button(panel, label="Run", pos=(280, 130))
        completeAnalysisButtonRun = wx.Button(panel, label="Complete analysis", pos=(550, 310))
        
        #CheckBoxes
        self.comparativeCheckBox = wx.CheckBox(panel, label="Comparative analysis", pos=(260, 110))
        self.autoCloseCheckBox = wx.CheckBox(panel, label="Close when finish", pos=(430, 315))
        #TextCtrls
        global fileInfo
        fileInfo = wx.TextCtrl(panel, -1, value="User: "+self.userName+"\r", pos = (10,10), size = (210, 400), style = wx.TE_MULTILINE)

        #Editables
        self.savePath.SetEditable(False)
        self.FileNameTextCtrl.SetEditable(False)
        self.saveButton.Disable()
        self.autoCloseCheckBox.Disable()
        mass3dButton.Disable()
        fileInfo.SetEditable(False)

        #Menubar
        fileButton = wx.Menu()
        helpButton = wx.Menu()
        aboutItem = helpButton.Append(wx.ID_ABOUT, 'About', 'About')
        openItem = fileButton.Append(wx.ID_OPEN, 'O&pen\tCtrl+O', 'Open chromatogram (.csv)')
        exitItem = fileButton.Append(wx.ID_EXIT, 'Q&uit\tCtrl+X', 'Exit the program')
        menuBar.Append(fileButton, 'File')
        menuBar.Append(helpButton, 'Help')

        self.SetMenuBar(menuBar)

        #Events
        self.Bind(wx.EVT_MENU, self.OnOpenFile, openItem)
        self.Bind(wx.EVT_BUTTON, self.OnOpenFile, self.openButton)
        self.Bind(wx.EVT_BUTTON, self.initialRunT, initialButtonRun)
        self.Bind(wx.EVT_MENU, self.Quit, exitItem)
        self.Bind(wx.EVT_CLOSE, self.Quit)
        self.Bind(wx.EVT_MENU, self.OnAboutBox, aboutItem)
        self.Bind(wx.EVT_BUTTON, self.cargarT, simpleAnalysisButtonRun)
        self.Bind(wx.EVT_BUTTON, self.completeAnalysisT, completeAnalysisButtonRun)
        self.Bind(wx.EVT_BUTTON, self.genCrom, cromatogramaButton)
        self.Bind(wx.EVT_BUTTON, self.massEspectrumForPeaksT, peaksMassSpectra)
        self.Bind(wx.EVT_BUTTON, self.mass3dEvt, mass3dButton)
        self.Bind(wx.EVT_BUTTON, self.OnSaveAs, self.saveButton)
        self.Bind(wx.EVT_BUTTON, self.OnOpenSpectralData, self.spectralDataButton)
        self.Bind(wx.EVT_BUTTON, self.cargarSpectralDataEvt, self.spectralDataButtonRun)
        self.Bind(wx.EVT_CHECKBOX, self.savePathHandler, self.automaticSaveCheckBox)
        self.SetTitle('Rutherford Chemistry: ' + self.userName)
        self.Show(True)

    def OnAboutBox(self, e):
        
        description = """Rutherford Chemistry is a Python-based GC-MS data analyzer for the Windows operating system. 
Features include a powerful built-in plot editor provided by Matplotlib; peak detection, and a Wiener filter both from SciPy; 
advanced compound search capabilities, file comparison and more.

It has been developed as a complementary software for the analytical instruments in Universidad de los Andes.
"""

        licence = """Rutherford Chemistry is free software; you can redistribute 
it and/or modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation.

Rutherford Chemistry is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
"""

        info = wx.AboutDialogInfo()

        info.SetIcon(wx.Icon('icon.ico', wx.BITMAP_TYPE_ICO))
        info.SetName('Rutherford Chemistry')
        info.SetVersion('1.0')
        info.SetDescription(description)
        info.SetCopyright('(C) 2015 Juan Barbosa')
        info.SetLicence(licence)
        
        wx.AboutBox(info)
           
    def savePathHandler(self, e):
        temp = self.automaticSaveCheckBox.GetValue()
        self.savePath.SetLabel("")
        self.savePath.SetEditable(temp)
        self.extention = ""
        if temp:
            self.saveButton.Enable()
            self.autoCloseCheckBox.Enable()
        else:
            self.saveButton.Disable()
            self.autoCloseCheckBox.Disable()

    def Quit(self, e):
        yesNoBox = wx.MessageDialog(None, 'Are you sure you want to continue?', 'Question', wx.YES_NO)
        yesNoAnswer = yesNoBox.ShowModal()
        yesNoBox.Destroy()
        if yesNoAnswer == wx.ID_YES:
            os._exit(1)
            
    def sysExit(self):
        sys.exit()
            
    def initialRunT(self, e):
        initialRunThread()
    
    def initialRun(self):
        wx.BeginBusyCursor()
        if len(self.filePaths) != 0:
            initialCromatogramasPos = len(self.cromatogramas)
            initialFilePathsPos = len(self.filePaths)
            for i in range(initialCromatogramasPos, initialFilePathsPos):
                temp = cromatogramaClass(numpy.genfromtxt(self.filePaths[i], delimiter=";"))
                self.cromatogramas.append(temp)
                fileInfo.AppendText("File "+str(i+1)+": "+ self.filePaths[i] + "\r")
            fileInfo.AppendText("--Chromatographic paths have been loaded-- \r")
        else:
            dlgM = wx.MessageDialog(self, "No file has been loaded", "Error", wx.OK)
            result = dlgM.ShowModal()
        wx.EndBusyCursor()

    def OnOpenFile(self, event):
        global fileInfo
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=self.currentDirectory, defaultFile="", wildcard=wildcard, style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            for path in paths:
                if path in self.filePaths:
                    dlgM = wx.MessageDialog(self, "File has been loaded already", "Error", wx.OK)
                    result = dlgM.ShowModal()
                else:
                    self.filePaths.append(path)
                    self.FileNameTextCtrl.AppendText(path+" ")
        dlg.Destroy()
        
    def OnSaveAs(self, event):
        saveFileDialog = wx.FileDialog(self, "Save directory", "", ".pdf", "Portable Document Format (*.pdf)|*.pdf|""Encapsulated Postscript (*.eps)|*.eps|""Joint Photographic Experts Group (*.jpge,*jpg)|*.jpg|""PGF code for LaTex (*.pgf)|*.pgf|""Raw RGBA bitmap (*.raw, *.rgba)|*.raw|""Scalable Vector Graphics (*.svg, *.svgz)|*.svg|""Tagged Image File Format (*tif, *.tiff)|*.tif", wx.FD_SAVE)
        if saveFileDialog.ShowModal() == wx.ID_CANCEL:
            return     # the user changed idea...
        else:
            temp = saveFileDialog.GetPath()[:-4]
            self.savePath.SetLabel(temp)
            self.extention = saveFileDialog.GetPath()[-4:]
            
    def genCrom(self, e):
        self.generarCromatograma()
    
    def generarCromatograma(self):
        global fileInfo
        wx.BeginBusyCursor()
        if len(self.cromatogramas) != 0 and self.comparativeCheckBox.GetValue() == False:
            for i in range(len(self.cromatogramas)):
                title = self.filePaths[i][:-4].rsplit('\\', 1)[-1]
                temp = self.savePath.GetLabel() + title +self.extention
                self.cromatogramas[i].generarCromatograma(self.userName, title, temp)
            fileInfo.AppendText("--Chromatographic data has been generated.-- \r")
        elif len(self.cromatogramas) != 0 and self.comparativeCheckBox.GetValue() == True:
            self.comparativeChromatogram()
            fileInfo.AppendText("Comparative chromatogram has been generated.\n")
            fileInfo.AppendText("--Chromatographic data has been generated.-- \r")
        else:
            dlgM = wx.MessageDialog(self, "No file has been loaded.", "Error", wx.OK)
            result = dlgM.ShowModal()
            fileInfo.AppendText("Error, no file has been loaded.\n")
        wx.EndBusyCursor()
        
    def cargarT(self, e):
        cargarThread()
        
    def cargar(self):
        global fileInfo
        wx.BeginBusyCursor()
        if len(self.cromatogramas) != 0 and self.comparativeCheckBox.GetValue() == False:
            for i in range(len(self.cromatogramas)):
                title = self.filePaths[i][:-4].rsplit('\\', 1)[-1]
                temp = self.savePath.GetLabel() + title +self.extention
                if self.cromatogramas[i].darChecker() == -1:
                    self.cromatogramas[i].inicializador()
            fileInfo.AppendText("--Chromatographic data has been loaded.-- \r")
        else:
            dlgM = wx.MessageDialog(self, "No file has been loaded.", "Error", wx.OK)
            result = dlgM.ShowModal()
            fileInfo.AppendText("Error, no file has been loaded.\n")
        wx.EndBusyCursor()
        
    def comparativeChromatogram(self):
        wx.BeginBusyCursor()
        temp = self.savePath.GetLabel() + "ComparativeAnalysis" + self.extention
        i = 0
        plt.figure(figsize=(10,5))
        for item in self.cromatogramas:
            labelText = self.filePaths[i][:-4].rsplit('\\', 1)[-1]
            plt.plot(item.darTiempo(), item.darIntensidad(), label=labelText)
            i += 1
        plt.legend()
        plt.xlabel("Tiempo (min)")
        plt.ylabel("Intensidad (u.a.)")
        plt.title("Comparative analysis")
        plt.suptitle(self.userName)
        if len(temp) == 19:
            plt.show()
        else:
            plt.savefig(temp)
            plt.close()
        wx.EndBusyCursor()
        
    def OnOpenSpectralData(self, event):
        global fileInfo
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=self.currentDirectory, defaultFile="", wildcard=wildcard2, style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            i = 1
            for path in paths:
                if path in self.spectralPaths:
                    dlgM = wx.MessageDialog(self, "File has been loaded already", "Error", wx.OK)
                    result = dlgM.ShowModal()
                else:
                    self.spectralPaths.append(path)
                    self.spectralDataTextCtrl.AppendText(path + " ")
                    fileInfo.AppendText("Spectral data: " + str(i) + " " + path +".\r")
            i += 1
            fileInfo.AppendText("--Spectral Data has been loaded--\r")
            
    def cargarSpectralDataEvt(self, e):
        self.cargarSpectralData()
            
    def cargarSpectralData(self):
        global fileInfo
        wx.BeginBusyCursor()
        if len(self.spectralPaths) != 0 and len(self.cromatogramas) != 0:
            for i in range(len(self.cromatogramas)):
                cromName = self.filePaths[i][:-4].rsplit('\\', 1)[-1]
                j = 0
                plt.figure(figsize=(10,5))
                automaticSaving = self.savePath.GetLabel() + cromName + "_SpectralData"
                for item in self.spectralPaths:
                    temp = numpy.genfromtxt(item)
                    fileName = item[:-4].rsplit('\\', 1)[-1]
                    tempMZ = []
                    tempInt = []
                    for fila in temp:
                        tempMZ.append(fila[0])
                        tempInt.append(fila[1])
                    self.cromatogramas[i].buscador(tempMZ, tempInt, 10, self.userName, fileName, j)
                    fileInfo.AppendText("File "+str(i+1)+" "+fileName+" has been processed \r")
                    j += 1
                    plt.xlabel("Tiempo (min)")
                    plt.ylabel("Intensidad (u.a.)")
                    plt.suptitle(self.userName)
                    plt.title(cromName)
                    plt.legend()
                if len(automaticSaving) == len(str(cromName) + "_SpectralData"):
                    plt.show()
                else:
                    plt.savefig(automaticSaving + self.extention)
                    plt.close()
            fileInfo.AppendText("--Spectral Data has been created-- \r")
        elif len(self.cromatogramas) == 0:
            dlgM = wx.MessageDialog(self, "No file has been loaded.", "Error", wx.OK)
            result = dlgM.ShowModal()
            fileInfo.AppendText("Error, no file has been loaded.\n")   
        else:
            fileInfo.AppendText("Spectral Data is not specified. \r")
        wx.EndBusyCursor()
        
    def massEspectrumForPeaksT(self, e):
        massEspectrumForPeaksThread()        
        
    def massEspectrumForPeaks(self):
        global fileInfo
        wx.BeginBusyCursor()
        if len(self.savePath.GetLabel()) != 0 and len(self.cromatogramas) != 0:
            i = 0
            for cromatograma in self.cromatogramas: 
                temp = self.savePath.GetLabel()+self.filePaths[i][:-4].rsplit('\\', 1)[-1]
                cromatograma.darEspectrosMasas(self.userName, temp, self.extention)
                fileInfo.AppendText(self.filePaths[i][:-4].rsplit('\\', 1)[-1] + " mass espectrum for peaks has been generated. \r")
                i += 1
        elif len(self.savePath.GetLabel()) == 0:
            dlgM = wx.MessageDialog(self, "No saving path is especified.", "Error", wx.OK)
            result = dlgM.ShowModal()
            fileInfo.AppendText("Error, no saving path is especified.\n")
        else:
            dlgM = wx.MessageDialog(self, "No file has been loaded.", "Error", wx.OK)
            result = dlgM.ShowModal()
            fileInfo.AppendText("Error, no file has been loaded.\n")
        wx.EndBusyCursor()
    
    def mass3dEvt(self, event):
        self.mass3d()
        
    def mass3d(self):
        global fileInfo
        fileInfo.AppendText("Still working on 3D mass spectrum" + "\r")
    
    def darAutoClose(self):
        return self.autoCloseCheckBox.GetValue()
        
    def completeAnalysis(self):
        self.initialRun()
        self.cargar()
        self.generarCromatograma()
        self.comparativeChromatogram()
        self.mass3d()
        self.massEspectrumForPeaks()
        self.cargarSpectralData()

    def completeAnalysisT(self, e):
        global fileInfo
        if len(self.filePaths) != 0 and len(self.savePath.GetLabel()) != 0:
            completeAnalysisThread()                            
        elif len(self.filePaths) == 0:
            dlgM = wx.MessageDialog(self, "No file has been loaded.", "Error", wx.OK)
            result = dlgM.ShowModal()
            fileInfo.AppendText("Error, no file has been loaded.\n")
        elif len(self.savePath.GetLabel()) == 0:
            dlgM = wx.MessageDialog(self, "No saving path is especified.", "Error", wx.OK)
            result = dlgM.ShowModal()
            fileInfo.AppendText("Error, no saving path is especified.\n")
"""
Clase cromatograma
"""
class cromatogramaClass():
    def __init__(self, matriz):
        global fileInfo
        self.checker = -1
        self.filas = []
        self.intensidad = []
        self.tiempo = []
        for i in range(1,len(matriz)):
            temp = filaClass(matriz[i], 2)
            self.filas.append(temp)
            if i%100 == 0:
                fileInfo.AppendText("Row " + str(i) + " has been created.\r")

    def calcularIntensidad(self):
        global fileInfo
        i = 0
        for item in self.filas:
            item.inicializador()
            self.intensidad.append(item.darSuma())
            if i%100 == 0:
                fileInfo.AppendText("Row " + str(i) + " intensity has been calculated.\r")
            i += 1

    def calcularTiempo(self):
        global fileInfo
        i = 0
        for item in self.filas:
            self.tiempo.append(item.darTiempo())
            if i%100 == 0:
                fileInfo.AppendText("Row " + str(i) + " time has been calculated.\r")
            i += 1
            
    def darTiempo(self):
        return self.tiempo

    def darIntensidad(self):
        return self.intensidad

    def inicializador(self):
        self.checker = 0
        self.calcularIntensidad()
        self.calcularTiempo()
        self.filterer()
        self.peakfinder()
    def darChecker(self):
        return self.checker
        
    def filterer(self):
        self.filtered = []
        global fileInfo
        temp = []
        temp.append(wiener(self.intensidad))
        for i in range (10):
            tempI = wiener(temp[i])
            temp.append(tempI)
            fileInfo.AppendText("Wiener filter pass " + str(i+1) + "\r")
        self.filtered = temp[-1]
    
    def peakfinder(self):
        global fileInfo
        nose = numpy.arange(0.1, 100)
        fileInfo.AppendText("Finding peaks... \r")
        self.peakind = find_peaks_cwt(self.filtered, nose, noise_perc=10)
        self.numberpeaks = len(self.peakind)
        self.labelText = str(self.numberpeaks) + " picos"
        self.peakX = []
        self.peakY = []
        fileInfo.AppendText("Appending peaks...\r")
        for peak in self.peakind:
            self.peakX.append(self.tiempo[peak])
            self.peakY.append(self.intensidad[peak]+10000000)

    def generarCromatograma(self, author, title, automaticSaving):
        global fileInfo
        fileInfo.AppendText("Generating Chromatogram...\r")
        plt.figure(figsize=(10,5))
        plt.plot(self.tiempo, self.intensidad, color="blue")
        plt.plot(self.peakX, self.peakY, "+", color="r", label=self.labelText)
        plt.xlabel("Tiempo (min)")
        plt.ylabel("Intensidad (u.a.)")
        plt.suptitle(author)
        plt.title(title)
        plt.legend()
        if len(automaticSaving) == len(title):
            fileInfo.AppendText(title + " chromatogram has been saved. \r")
            plt.show()
        else:
            plt.savefig(automaticSaving)
            plt.close()
            fileInfo.AppendText(title + " chromatogram has been saved.\r")

    def darEspectrosMasas(self, author, automaticSaving, extention):
        global fileInfo
        for peak in self.peakind:
            tempStr = str(self.tiempo[peak])
            if str(self.tiempo[peak]) >= 5:
                tempStr = str(self.tiempo[peak])[:5]
            temp = automaticSaving + "_Tiempo_" + tempStr + extention
            self.filas[peak].espectroMasas(author, temp)
            fileInfo.AppendText("Mass spectrum for retention time " + tempStr + " min has been saved. \r")

    def buscador(self, mz, intensidad, acotado, author, title, pos):
        aciertos = []
        y = []
        i = 0
        for fila in self.filas:
            aciertos.append(fila.comparador(mz, intensidad, acotado))
        for acierto in aciertos:
            if acierto != 0:
                y.append(self.intensidad[i])
            else:
                y.append(0)
            i += 1
        if pos == 0:
            plt.plot(self.tiempo, self.intensidad, linewidth = 0.1, color='0.8')
        base_line, = plt.plot([], [], linewidth=10, label=title)
        plt.fill_between(self.tiempo, 0, y, facecolor=base_line.get_color(), edgecolor="none")

class filaClass(cromatogramaClass):
    def __init__(self, lista, beginColumn):
        self.filaCompleta = []
        self.workingFila = []
        self.suma = 0
        self.maximo = 0
        self.posicion = -1
        self.relacionMZ = []
        self.tiempo = 0
        self.beginColumn = beginColumn
        for item in lista:
            self.filaCompleta.append(item)

    def filaBeginColumn(self):
        for i in range(self.beginColumn, len(self.filaCompleta)):
            self.workingFila.append(self.filaCompleta[i])

    def calcularIntensidadRelativa(self):
        self.intensidadRelativa = []
        maxIntensity = max(self.workingFila)
        for item in self.workingFila:
            if maxIntensity != 0:
                temp = item/maxIntensity * 100
                self.intensidadRelativa.append(temp)
            else:
                self.intensidadRelativa.append(0)

    def darWorkingFila(self):
        return self.workingFila

    def sumador(self):
        for item in self.workingFila:
            self.suma += item

    def darSuma(self):
        return self.suma

    def valorMaximo(self):
        self.maximo = max(self.workingFila)

    def darValorMaximo(self):
        return self.maximo

    def posicionMaximo(self):
        self.posicion = self.workingFila.index(self.maximo)

    def darPosicionMaximo(self):
        return self.posicion

    def relacion(self):
        for i in range(len(self.workingFila)):
            self.relacionMZ.append(i+1)

    def calcularRelacionMZshorther(self):
        temp = []
        for i in range(len(self.relacionMZ)):
            if self.workingFila[i] != 0:
                temp.append(i)
        self.maxI = max(temp)

    def espectroMasas(self, author, automaticSaving):
        self.calcularRelacionMZshorther()
        temp = self.maxI + 1
        plt.figure(figsize=(10,5))
        plt.bar(self.relacionMZ[:temp], self.intensidadRelativa[:temp], width=0.2, color="black")
        plt.title(automaticSaving[:-4].rsplit('_', 1)[-1] + " min")
        plt.xlabel("m/z")
        plt.ylabel("Intensidad relativa")
        plt.suptitle(author)
        plt.savefig(automaticSaving)
        plt.close()

    def comparador(self, mz, intensity, acotado):
        aciertos = 0
        total = len(mz)
        maxintR = mz[intensity.index(max(intensity))]
        maxintI = self.relacionMZ[self.intensidadRelativa.index(max(self.intensidadRelativa))]
        if maxintR != maxintI:
            return 0
        else:
            for i in range (len(mz)):
                pos = self.relacionMZ.index(mz[i])
                expIntensity = self.intensidadRelativa[pos]
                patterIntensity = intensity[i]
                minPossible = patterIntensity - acotado
                maxPossible = patterIntensity + acotado
                if expIntensity > minPossible and expIntensity <= maxPossible:
                    aciertos += 1.0
                    return (aciertos/total)*100

    def calcularTiempo(self):
        self.tiempo = self.filaCompleta[1]

    def darTiempo(self):
        return self.tiempo

    def inicializador(self):
        self.filaBeginColumn()
        self.sumador()
        self.valorMaximo()
        self.posicionMaximo()
        self.relacion()
        self.calcularTiempo()
        self.calcularIntensidadRelativa()

def main():
    app = wx.App()
    global mainWindow
    mainWindow = windowClass(None)
    app.MainLoop()

main()