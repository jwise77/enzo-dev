# Root wrapper Makefile for CMake transition

BUILD_DIR ?= build

.PHONY: all build clean install enzo inits enzohop ring anyl P-GroupFinder

all: build

build:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --parallel

clean:
	@rm -rf $(BUILD_DIR)

enzo:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target enzo --parallel

inits:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target inits --parallel

enzohop:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target enzohop --parallel

ring:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target ring --parallel

anyl:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target anyl --parallel

P-GroupFinder:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target P-GroupFinder --parallel
