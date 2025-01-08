PROJECT_NAME = main

CXXFLAGS += -DDEBUG

INPUT_FILE = input.txt
OUTPUT_FILE = output.txt

NPROCS ?= $(shell nproc)

BUILD_DIR = build
RESULTS_DIR = results

all: build test run

build:
	@echo "==> Configuring the project..."
	@cmake -B$(BUILD_DIR) -H.
	@echo "==> Building the project..."
	@cmake --build $(BUILD_DIR) -j $(NPROCS)

test: build
	@echo "==> Running tests..."
	@cd $(BUILD_DIR) && ctest --verbose

run: build
	@echo "==> Running ${PROJECT_NAME} with arguments: $(ARGS)"
	@${BUILD_DIR}/${PROJECT_NAME} $(ARGS)

run-all: build
	@echo "==> Running ${PROJECT_NAME} for different amount of elements"
	@python3 solve.py

clean:
	@echo "==> Cleaning up..."
	@make -C ${BUILD_DIR} clean

install:
	sudo apt-get update
	sudo apt-get install -y cmake clang libgtest-dev

update:
	git submodule update --remote --merge --recursive

.PHONY: all build test run run-all clean install update
