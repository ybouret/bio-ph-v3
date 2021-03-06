CREATE := ./create.sh

all:

clean:
	@echo "-- removing temporary files" && rm -f *.bin *.dat *.vtk *.csv
	@echo "-- removing junk..." && find . -name '*~' | xargs rm -f


veryclean: clean
	@echo "-- removing out of sources builds" && cd forge && touch targets && ( cat targets | xargs rm -rf ) && rm -f targets
	@echo "-- removing local binaries" && rm -Rf ./bin

opt:
	@bash ./build.sh clang Release

gnu:
	@bash $(CREATE) src gnu ${BUILD_TYPE}

intel:
	@bash $(CREATE) src intel ${BUILD_TYPE}

clang:
	@bash $(CREATE) src clang ${BUILD_TYPE}

xcode:
	@bash $(CREATE) src xcode ${BUILD_TYPE}
	
vs9:
	@bash $(CREATE) src vs9 ${BUILD_TYPE}

vs10:
	@bash $(CREATE) src vs10 ${BUILD_TYPE}

