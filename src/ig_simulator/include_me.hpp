#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "stdlib.h"
#include "stdio.h"
#include <memory>
#include "../utils/sequence_tools.hpp"

size_t RandomInt(size_t max) {
    return rand() % max;
}

size_t RandomIndex(size_t to, size_t from = 0) {
    return rand() % (to - from) + from;
}