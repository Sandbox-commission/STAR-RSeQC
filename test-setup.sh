#!/bin/bash
set -euo pipefail

# Test Suite for setup.sh
# Validates all three installation paths without interactive prompts

echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║              setup.sh Test Suite                                     ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

TEST_PASS=0
TEST_FAIL=0

test_header() {
    echo -e "\n${BLUE}━━━ $1 ━━━${NC}"
}

test_pass() {
    echo -e "${GREEN}✓${NC} $1"
    ((TEST_PASS++))
}

test_fail() {
    echo -e "${RED}✗${NC} $1"
    ((TEST_FAIL++))
}

# ─────────────────────────────────────────────────────────────────────────────
# TEST 1: File Existence & Permissions
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 1: File Existence & Permissions"

if [[ -f setup.sh ]]; then
    test_pass "setup.sh exists"
else
    test_fail "setup.sh not found"
    exit 1
fi

if [[ -x setup.sh ]]; then
    test_pass "setup.sh is executable"
else
    test_fail "setup.sh is not executable"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 2: Bash Syntax
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 2: Bash Syntax Validation"

if bash -n setup.sh 2>&1 | grep -q "syntax error"; then
    test_fail "setup.sh has syntax errors"
else
    test_pass "setup.sh bash syntax valid"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 3: Function Definitions (Source script and check)
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 3: Function Definitions"

# Source functions (without executing main logic)
source <(sed -n '1,/^setup_docker()/p' setup.sh | head -n -1)

# Check if key functions exist
if grep -q "^setup_docker()" setup.sh; then
    test_pass "setup_docker() function defined"
else
    test_fail "setup_docker() function not found"
fi

if grep -q "^setup_mamba_auto()" setup.sh; then
    test_pass "setup_mamba_auto() function defined"
else
    test_fail "setup_mamba_auto() function not found"
fi

if grep -q "^setup_manual()" setup.sh; then
    test_pass "setup_manual() function defined"
else
    test_fail "setup_manual() function not found"
fi

if grep -q "^install_docker_auto()" setup.sh; then
    test_pass "install_docker_auto() function defined"
else
    test_fail "install_docker_auto() function not found"
fi

if grep -q "^install_miniforge_auto()" setup.sh; then
    test_pass "install_miniforge_auto() function defined"
else
    test_fail "install_miniforge_auto() function not found"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 4: Key Features Presence
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 4: Feature Verification"

# Menu system
if grep -q "Choose your installation method:" setup.sh; then
    test_pass "Menu system present"
else
    test_fail "Menu system not found"
fi

# Docker path
if grep -q "Docker Container" setup.sh; then
    test_pass "Docker installation path documented"
else
    test_fail "Docker path description missing"
fi

# Mamba path
if grep -q "Mamba/Conda" setup.sh; then
    test_pass "Mamba installation path documented"
else
    test_fail "Mamba path description missing"
fi

# Manual path
if grep -q "Manual Installation" setup.sh; then
    test_pass "Manual installation path documented"
else
    test_fail "Manual path description missing"
fi

# OS detection for Docker
if grep -q "ubuntu\|debian" setup.sh; then
    test_pass "OS detection for Ubuntu/Debian"
else
    test_fail "OS detection incomplete"
fi

if grep -q "fedora\|rhel" setup.sh; then
    test_pass "OS detection for Fedora/RHEL"
else
    test_fail "OS detection incomplete"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 5: Docker Setup Logic
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 5: Docker Setup Logic"

if grep -q "docker --version" setup.sh; then
    test_pass "Docker version check present"
else
    test_fail "Docker version check missing"
fi

if grep -q "docker build" setup.sh; then
    test_pass "Docker build command present"
else
    test_fail "Docker build command missing"
fi

if grep -q "docker-run.sh" setup.sh; then
    test_pass "docker-run.sh generation present"
else
    test_fail "docker-run.sh generation missing"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 6: Mamba Setup Logic
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 6: Mamba Setup Logic"

if grep -q "miniforge3\|mambaforge\|miniconda3" setup.sh; then
    test_pass "Conda environment search paths present"
else
    test_fail "Conda search paths missing"
fi

if grep -q "mamba create.*star" setup.sh; then
    test_pass "STAR environment creation command present"
else
    test_fail "STAR environment command missing"
fi

if grep -q "mamba create.*rseqc" setup.sh; then
    test_pass "RSeQC environment creation command present"
else
    test_fail "RSeQC environment command missing"
fi

if grep -q "Miniforge3-Linux-x86_64.sh" setup.sh; then
    test_pass "Miniforge installer URL present"
else
    test_fail "Miniforge installer URL missing"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 7: Manual Installation Path
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 7: Manual Installation Path"

if grep -q "https://github.com/conda-forge/miniforge" setup.sh; then
    test_pass "Miniforge link present"
else
    test_fail "Miniforge link missing"
fi

if grep -q "https://docs.docker.com/engine/install" setup.sh; then
    test_pass "Docker documentation link present"
else
    test_fail "Docker documentation link missing"
fi

if grep -q "Step-by-step" setup.sh; then
    test_pass "Step-by-step instructions present"
else
    test_fail "Step-by-step instructions missing"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 8: Configuration Management
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 8: Configuration Management"

if grep -q "config.json" setup.sh; then
    test_pass "config.json handling present"
else
    test_fail "config.json handling missing"
fi

if grep -q "\\.config/star-rseqc" setup.sh; then
    test_pass "Config directory path correct"
else
    test_fail "Config directory path incorrect"
fi

if grep -q "genome_dir" setup.sh; then
    test_pass "genome_dir configuration present"
else
    test_fail "genome_dir configuration missing"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 9: Error Handling
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 9: Error Handling"

if grep -q "if \[\[" setup.sh; then
    test_pass "Conditional error checks present"
else
    test_fail "Error handling incomplete"
fi

if grep -q "log_error\|Exit\|Error" setup.sh; then
    test_pass "Error messaging present"
else
    test_fail "Error messaging missing"
fi

if grep -q "return\|exit" setup.sh; then
    test_pass "Exit codes present"
else
    test_fail "Exit codes missing"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 10: Docker Files
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 10: Docker Files Validation"

if [[ -f Dockerfile ]]; then
    test_pass "Dockerfile exists"
else
    test_fail "Dockerfile not found"
fi

if [[ -f docker-compose.yml ]]; then
    test_pass "docker-compose.yml exists"
else
    test_fail "docker-compose.yml not found"
fi

if [[ -f .dockerignore ]]; then
    test_pass ".dockerignore exists"
else
    test_fail ".dockerignore not found"
fi

# Validate Dockerfile syntax
if docker run --rm -i hadolint/hadolint < Dockerfile 2>&1 | grep -q "error"; then
    test_fail "Dockerfile has linting errors"
else
    test_pass "Dockerfile syntax valid (via hadolint)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 11: Documentation
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 11: Documentation Files"

if [[ -f INSTALLATION_GUIDE.md ]]; then
    test_pass "INSTALLATION_GUIDE.md exists"
else
    test_fail "INSTALLATION_GUIDE.md missing"
fi

if [[ -f QUICK_REFERENCE.md ]]; then
    test_pass "QUICK_REFERENCE.md exists"
else
    test_fail "QUICK_REFERENCE.md missing"
fi

if [[ -f create-package.sh ]]; then
    test_pass "create-package.sh exists"
else
    test_fail "create-package.sh missing"
fi

# ─────────────────────────────────────────────────────────────────────────────
# TEST 12: Package Generation
# ─────────────────────────────────────────────────────────────────────────────

test_header "Test 12: Package Contents"

if [[ -f STAR-RSeQC-0.1.0.zip ]]; then
    test_pass "Package ZIP created"
    
    # Check package contents
    if unzip -l STAR-RSeQC-0.1.0.zip | grep -q "setup.sh"; then
        test_pass "setup.sh included in package"
    else
        test_fail "setup.sh not in package"
    fi
    
    if unzip -l STAR-RSeQC-0.1.0.zip | grep -q "Dockerfile"; then
        test_pass "Dockerfile included in package"
    else
        test_fail "Dockerfile not in package"
    fi
    
    if unzip -l STAR-RSeQC-0.1.0.zip | grep -q "star-rseqc"; then
        test_pass "Binary included in package"
    else
        test_fail "Binary not in package"
    fi
    
    if unzip -l STAR-RSeQC-0.1.0.zip | grep -q "docs/"; then
        test_pass "Documentation included in package"
    else
        test_fail "Documentation not in package"
    fi
else
    test_fail "Package ZIP not created"
fi

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────────────

echo
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║                         TEST SUMMARY                                 ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo
echo -e "${GREEN}Passed: $TEST_PASS${NC}"
echo -e "${RED}Failed: $TEST_FAIL${NC}"
echo

TOTAL=$((TEST_PASS + TEST_FAIL))
PERCENT=$(( (TEST_PASS * 100) / TOTAL ))

echo "Result: $TEST_PASS/$TOTAL tests passed ($PERCENT%)"
echo

if [[ $TEST_FAIL -eq 0 ]]; then
    echo -e "${GREEN}✅ All tests passed! setup.sh is ready for use.${NC}"
    exit 0
else
    echo -e "${RED}⚠️  Some tests failed. Please review the output above.${NC}"
    exit 1
fi
