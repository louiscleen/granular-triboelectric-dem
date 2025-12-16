# --- Compilateur et répertoires ---
CXX := g++
SRC_DIR := src
BUILD_DIR := build
EXTERNAL_DIR := external
DATA_DIR := data
TARGET := dem

# --- Dossier d'installation des dépendances externes ---
VORO_DIR    := ../voro++-0.4.6/src
EIGEN_DIR   := $(EXTERNAL_DIR)/eigen-5.0.0
CXXOPTS_DIR := $(EXTERNAL_DIR)/cxxopts
TOML_DIR    := $(EXTERNAL_DIR)/toml




# --- Flags ---
CXXFLAGS := -std=c++20 -O3 -Wall -Wextra -Wshadow -flto \
			-Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR) -MMD -MP
LDFLAGS := -flto -L$(VORO_DIR) -lvoro++

CXXFLAGS_DEBUG := -std=c++20 -g -O0 -Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR)
LDFLAGS_DEBUG := -L$(VORO_DIR) -lvoro++

CXX_MPI := g++
CXXFLAGS_MPI := -std=c++20 -O3 -Wall -Wextra -Wshadow \
			-Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR) -MMD -MP
LDFLAGS_MPI := -L$(VORO_DIR) -lvoro++

CXXFLAGS_QUICK := -std=c++20 -O1 -Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -MMD -MP
LDFLAGS_QUICK :=

CXXFLAGS_FAST := -std=c++20 -O3 -march=native -mtune=native -fno-math-errno -fno-trapping-math -ffp-contract=fast -DNDEBUG \
			-Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR)
LDFLAGS_FAST := -L$(VORO_DIR) -lvoro++

CXXFLAGS_MPI_FAST := -std=c++20 -O3 -march=native -mtune=native -fno-math-errno -fno-trapping-math -ffp-contract=fast -DNDEBUG \
			-Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR)
LDFLAGS_MPI_FAST := -L$(VORO_DIR) -lvoro++

CXXFLAGS_LM4 := -std=c++20 -O3 -march=znver4 -mtune=znver4 -fno-math-errno -fno-trapping-math -ffp-contract=fast -DNDEBUG \
			-Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR)
LDFLAGS_LM4 := -L$(VORO_DIR) -lvoro++

CXXFLAGS_NIC5 := -std=c++20 -O3 -march=znver2 -mtune=znver2 -fno-math-errno -fno-trapping-math -ffp-contract=fast -DNDEBUG \
			-Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR)
LDFLAGS_NIC5 := -L$(VORO_DIR) -lvoro++

CXXFLAGS_D2 := -std=c++20 -O3 -march=skylake -mtune=skylake -fno-math-errno -fno-trapping-math -ffp-contract=fast -DNDEBUG \
			-Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR)
LDFLAGS_D2 := -L$(VORO_DIR) -lvoro++

CXXFLAGS_H2 := -std=c++20 -O3 -march=znver1 -mtune=znver1 -fno-math-errno -fno-trapping-math -ffp-contract=fast -DNDEBUG \
			-Iinclude -I$(EIGEN_DIR) -I$(CXXOPTS_DIR) -I$(TOML_DIR) -I$(VORO_DIR) -isystem $(VORO_DIR)
LDFLAGS_H2 := -L$(VORO_DIR) -lvoro++



# --- Fichiers sources et objets ---
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS_RELEASE := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/release/%.o)
OBJS_DEBUG := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/debug/%.o)
OBJS_MPI := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/MPI/%.o)
OBJS_QUICK := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/quick/%.o)
OBJS_FAST := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/fast/%.o)
OBJS_MPI_FAST := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/MPI_fast/%.o)
OBJS_LM4 := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/LM4/%.o)
OBJS_NIC5 := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/NIC5/%.o)
OBJS_D2 := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/D2/%.o)
OBJS_H2 := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/H2/%.o)

DEPS_RELEASE := $(OBJS_RELEASE:.o=.d)
DEPS_DEBUG := $(OBJS_DEBUG:.o=.d)
DEPS_MPI := $(OBJS_MPI:.o=.d)
DEPS_QUICK := $(OBJS_QUICK:.o=.d)
DEPS_FAST := $(OBJS_FAST:.o=.d)
DEPS_MPI_FAST := $(OBJS_MPI_FAST:.o=.d)

# --- MPI ---
# Détection Windows MS-MPI
ifeq ($(OS),Windows_NT)
	MSMPI_INC  := $(shell cygpath -u "$$MSMPI_INC")
	MSMPI_LIB64 := $(shell cygpath -u "$$MSMPI_LIB64")
	CXXFLAGS_MPI += -I"${MSMPI_INC}"
	LDFLAGS_MPI  += -L"${MSMPI_LIB64}" -lmsmpi
	CXXFLAGS_MPI_FAST += -I"${MSMPI_INC}"
	LDFLAGS_MPI_FAST += -L"${MSMPI_LIB64}" -lmsmpi
	CXXFLAGS_DEBUG += -I"${MSMPI_INC}"
	LDFLAGS_DEBUG += -L"${MSMPI_LIB64}" -lmsmpi

else
	# Linux / WSL / Mac → utiliser mpic++
	CXX_MPI := mpicxx
endif



# --- Règle par défaut : release ---
all: release

# --- Mode Release ---
release: $(TARGET)

$(TARGET): $(OBJS_RELEASE)
	$(CXX) $^ -o $@ $(LDFLAGS)  

$(BUILD_DIR)/release/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/release
	$(CXX) $(CXXFLAGS) -c $< -o $@

# --- Mode Debug ---
debug: $(TARGET)_debug

$(TARGET)_debug: $(OBJS_DEBUG)
	$(CXX) $^ -o $@ $(LDFLAGS_DEBUG)

$(BUILD_DIR)/debug/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/debug
	$(CXX) $(CXXFLAGS_DEBUG) -c $< -o $@

# --- Mode MPI ---
MPI: $(TARGET)_MPI

$(TARGET)_MPI: $(OBJS_MPI)
	$(CXX_MPI) $^ -o $@ $(LDFLAGS_MPI)  

$(BUILD_DIR)/MPI/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/MPI
	$(CXX_MPI) $(CXXFLAGS_MPI) -c $< -o $@

# --- Mode quick (pour tester rapidment) ---
quick: $(TARGET)_quick

$(TARGET)_quick: $(OBJS_QUICK)
	$(CXX) $^ -o $@ $(LDFLAGS_QUICK)  

$(BUILD_DIR)/quick/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/quick
	$(CXX) $(CXXFLAGS_QUICK) -c $< -o $@

# --- Mode fast ---
fast: $(TARGET)_fast

$(TARGET)_fast: $(OBJS_FAST)
	$(CXX) $^ -o $@ $(LDFLAGS_FAST)

$(BUILD_DIR)/fast/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/fast
	$(CXX) $(CXXFLAGS_FAST) -c $< -o $@

# --- Mode MPI_fast ---
MPI_fast: $(TARGET)_MPI_fast

$(TARGET)_MPI_fast: $(OBJS_MPI_FAST)
	$(CXX_MPI) $^ -o $@ $(LDFLAGS_MPI_FAST)

$(BUILD_DIR)/MPI_fast/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/MPI_fast
	$(CXX_MPI) $(CXXFLAGS_MPI_FAST) -c $< -o $@

# --- Mode LM4 ---
LM4: $(TARGET)_LM4

$(TARGET)_LM4: $(OBJS_LM4)
	$(CXX_MPI) $^ -o $@ $(LDFLAGS_LM4)

$(BUILD_DIR)/LM4/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/LM4
	$(CXX_MPI) $(CXXFLAGS_LM4) -c $< -o $@

# --- Mode NIC5 ---
NIC5: $(TARGET)_NIC5

$(TARGET)_NIC5: $(OBJS_NIC5)
	$(CXX_MPI) $^ -o $@ $(LDFLAGS_NIC5)

$(BUILD_DIR)/NIC5/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/NIC5
	$(CXX_MPI) $(CXXFLAGS_NIC5) -c $< -o $@

# --- Mode Dragon2 ---
D2: $(TARGET)_D2

$(TARGET)_D2: $(OBJS_D2)
	$(CXX_MPI) $^ -o $@ $(LDFLAGS_D2)
$(BUILD_DIR)/D2/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/D2
	$(CXX_MPI) $(CXXFLAGS_D2) -c $< -o $@

# --- Mode Hercule2 ---
H2: $(TARGET)_H2

$(TARGET)_H2: $(OBJS_H2)
	$(CXX_MPI) $^ -o $@ $(LDFLAGS_H2)
$(BUILD_DIR)/H2/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/H2
	$(CXX_MPI) $(CXXFLAGS_H2) -c $< -o $@

# --- Création des répertoires de build ---
$(BUILD_DIR)/release:
	mkdir -p $(BUILD_DIR)/release

$(BUILD_DIR)/debug:
	mkdir -p $(BUILD_DIR)/debug

$(BUILD_DIR)/MPI:
	mkdir -p $(BUILD_DIR)/MPI

$(BUILD_DIR)/quick:
	mkdir -p $(BUILD_DIR)/quick

$(BUILD_DIR)/fast:
	mkdir -p $(BUILD_DIR)/fast

$(BUILD_DIR)/MPI_fast:
	mkdir -p $(BUILD_DIR)/MPI_fast

$(BUILD_DIR)/LM4:
	mkdir -p $(BUILD_DIR)/LM4

$(BUILD_DIR)/NIC5:
	mkdir -p $(BUILD_DIR)/NIC5

$(BUILD_DIR)/D2:
	mkdir -p $(BUILD_DIR)/D2

$(BUILD_DIR)/H2:
	mkdir -p $(BUILD_DIR)/H2

# --- Nettoyage ---
clean:
	rm -rf $(BUILD_DIR) $(TARGET) $(TARGET)_debug $(TARGET)_MPI $(TARGET)_quick $(TARGET)_fast $(TARGET)_MPI_fast $(TARGET)_LM4 $(TARGET)_NIC5 $(TARGET)_D2 $(TARGET)_H2

# -- Inclure les dépendances automatiquement ---
-include $(DEPS_RELEASE)
-include $(DEPS_DEBUG)
-include $(DEPS_MPI)
-include $(DEPS_QUICK)
-include $(DEPS_FAST)

# --- Suppression des données ---
remove:
	rm -rf $(DATA_DIR)

.PHONY: all clean release debug mpi quick fast remove