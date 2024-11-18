PROJECT_NAME = main

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
	@echo "==> Running ${PROJECT_NAME}"
	@cd ${BUILD_DIR} && ./${PROJECT_NAME}

clean:
	@echo "==> Cleaning up..."
	@rm -rf $(BUILD_DIR)

rebuild: clean build

install:
	sudo apt-get update
	sudo apt-get install -y cmake clang libgtest-dev

.PHONY: all build test clean rebuild install
