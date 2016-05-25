CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

%.o : %.cc
	$(CXX) -c $<

all: polRapPtPlot polGenRecFitPlot EvaluateSyst PlotFinalResults MattsResults AverageSystematics ChangeTGraph AlterPPD
#all: polRapPtPlot polGenRecFitPlot PlotFinalResults

polRapPtPlot: polRapPtPlot.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

polGenRecFitPlot: polGenRecFitPlot.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

EvaluateSyst: EvaluateSyst.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

PlotFinalResults: PlotFinalResults.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

MattsResults: MattsResults.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

AverageSystematics: AverageSystematics.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

ChangeTGraph: ChangeTGraph.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

AlterPPD: AlterPPD.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit
									
clean: 
	rm polRapPtPlot polGenRecFitPlot EvaluateSyst PlotFinalResults MattsResults AverageSystematics ChangeTGraph AlterPPD*.o
