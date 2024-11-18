PROJECT_NAME = matrix

INPUT_FILE = input.txt
OUTPUT_FILE = output.txt

NPROCS ?= $(shell nproc)

BUILD_DIR = build

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
	@cd ${BUILD_DIR} && ./${PROJECT_NAME} < ../${INPUT_FILE} > ../${OUTPUT_FILE}

clean:
	@echo "==> Cleaning up..."
	@rm -rf $(BUILD_DIR)

rebuild: clean build

install:
	sudo apt-get update
	sudo apt-get install -y cmake clang libgtest-dev

.PHONY: all build test clean rebuild install
