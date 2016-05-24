##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=Hummingbird
ConfigurationName      :=Debug
WorkspacePath          := "/home/kp/Uni/NumerischeStroemungssimulation/codeliteWorkspace"
ProjectPath            := "/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird"
IntermediateDirectory  :=Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=KP
Date                   :=24/05/16
CodeLitePath           :="/home/kp/.codelite"
LinkerName             :=/usr/bin/g++-4.9
SharedObjectLinkerName :=/usr/bin/g++-4.9 -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird/Debug/Hummingbird.exe
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="Hummingbird.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++-4.9
CC       := /usr/bin/gcc-4.9
CXXFLAGS :=  -g $(Preprocessors)
CFLAGS   :=  -g $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/src_data.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_input.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_main.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_metrics.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_output.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_setup.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_solve.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d Debug || $(MakeDirCommand) Debug


$(IntermediateDirectory)/.d:
	@test -d Debug || $(MakeDirCommand) Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/src_data.cpp$(ObjectSuffix): src/data.cpp $(IntermediateDirectory)/src_data.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird/src/data.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_data.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_data.cpp$(DependSuffix): src/data.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_data.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_data.cpp$(DependSuffix) -MM "src/data.cpp"

$(IntermediateDirectory)/src_data.cpp$(PreprocessSuffix): src/data.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_data.cpp$(PreprocessSuffix) "src/data.cpp"

$(IntermediateDirectory)/src_input.cpp$(ObjectSuffix): src/input.cpp $(IntermediateDirectory)/src_input.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird/src/input.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_input.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_input.cpp$(DependSuffix): src/input.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_input.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_input.cpp$(DependSuffix) -MM "src/input.cpp"

$(IntermediateDirectory)/src_input.cpp$(PreprocessSuffix): src/input.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_input.cpp$(PreprocessSuffix) "src/input.cpp"

$(IntermediateDirectory)/src_main.cpp$(ObjectSuffix): src/main.cpp $(IntermediateDirectory)/src_main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird/src/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_main.cpp$(DependSuffix): src/main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_main.cpp$(DependSuffix) -MM "src/main.cpp"

$(IntermediateDirectory)/src_main.cpp$(PreprocessSuffix): src/main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_main.cpp$(PreprocessSuffix) "src/main.cpp"

$(IntermediateDirectory)/src_metrics.cpp$(ObjectSuffix): src/metrics.cpp $(IntermediateDirectory)/src_metrics.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird/src/metrics.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_metrics.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_metrics.cpp$(DependSuffix): src/metrics.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_metrics.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_metrics.cpp$(DependSuffix) -MM "src/metrics.cpp"

$(IntermediateDirectory)/src_metrics.cpp$(PreprocessSuffix): src/metrics.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_metrics.cpp$(PreprocessSuffix) "src/metrics.cpp"

$(IntermediateDirectory)/src_output.cpp$(ObjectSuffix): src/output.cpp $(IntermediateDirectory)/src_output.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird/src/output.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_output.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_output.cpp$(DependSuffix): src/output.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_output.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_output.cpp$(DependSuffix) -MM "src/output.cpp"

$(IntermediateDirectory)/src_output.cpp$(PreprocessSuffix): src/output.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_output.cpp$(PreprocessSuffix) "src/output.cpp"

$(IntermediateDirectory)/src_setup.cpp$(ObjectSuffix): src/setup.cpp $(IntermediateDirectory)/src_setup.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird/src/setup.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_setup.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_setup.cpp$(DependSuffix): src/setup.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_setup.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_setup.cpp$(DependSuffix) -MM "src/setup.cpp"

$(IntermediateDirectory)/src_setup.cpp$(PreprocessSuffix): src/setup.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_setup.cpp$(PreprocessSuffix) "src/setup.cpp"

$(IntermediateDirectory)/src_solve.cpp$(ObjectSuffix): src/solve.cpp $(IntermediateDirectory)/src_solve.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/kp/Uni/NumerischeStroemungssimulation/Hummingbird/src/solve.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_solve.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_solve.cpp$(DependSuffix): src/solve.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_solve.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_solve.cpp$(DependSuffix) -MM "src/solve.cpp"

$(IntermediateDirectory)/src_solve.cpp$(PreprocessSuffix): src/solve.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_solve.cpp$(PreprocessSuffix) "src/solve.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r Debug/


