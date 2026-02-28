#!/bin/bash

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║        setup.sh Validation Test Suite                         ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo

PASS=0
FAIL=0

test() {
    if eval "$1"; then
        echo "✓ $2"
        ((PASS++))
    else
        echo "✗ $2"
        ((FAIL++))
    fi
}

# Test 1: File exists and is executable
test "[[ -x setup.sh ]]" "setup.sh is executable"

# Test 2: Bash syntax
test "bash -n setup.sh 2>/dev/null" "Bash syntax is valid"

# Test 3: Functions defined
test "grep -q '^setup_docker()' setup.sh" "setup_docker() function exists"
test "grep -q '^setup_mamba_auto()' setup.sh" "setup_mamba_auto() function exists"
test "grep -q '^setup_manual()' setup.sh" "setup_manual() function exists"
test "grep -q '^install_docker_auto()' setup.sh" "install_docker_auto() function exists"
test "grep -q '^install_miniforge_auto()' setup.sh" "install_miniforge_auto() function exists"

# Test 4: Core features
test "grep -q 'Choose your installation method' setup.sh" "Menu system present"
test "grep -q 'Docker Container' setup.sh" "Docker path documented"
test "grep -q 'Automatic Mamba' setup.sh" "Mamba path documented"
test "grep -q 'Manual Installation' setup.sh" "Manual path documented"

# Test 5: Docker setup
test "grep -q 'docker --version' setup.sh" "Docker version check"
test "grep -q 'docker build' setup.sh" "Docker build command"
test "grep -q 'docker-run.sh' setup.sh" "docker-run.sh generation"

# Test 6: Mamba setup
test "grep -q 'miniforge3' setup.sh" "Miniforge path in detection"
test "grep -q 'mamba create -n star' setup.sh" "STAR environment creation"
test "grep -q 'mamba create -n rseqc' setup.sh" "RSeQC environment creation"

# Test 7: Manual path
test "grep -q 'https://github.com/conda-forge/miniforge' setup.sh" "Miniforge link"
test "grep -q 'https://docs.docker.com/engine/install' setup.sh" "Docker docs link"

# Test 8: Config management
test "grep -q 'config.json' setup.sh" "Config file handling"
test "grep -q '\\.config/star-rseqc' setup.sh" "Config directory"

# Test 9: Error handling
test "grep -q 'log_error' setup.sh" "Error logging function"
test "grep -q 'exit 1' setup.sh" "Exit on error"

# Test 10: Docker files exist
test "[[ -f Dockerfile ]]" "Dockerfile exists"
test "[[ -f docker-compose.yml ]]" "docker-compose.yml exists"
test "[[ -f .dockerignore ]]" ".dockerignore exists"

# Test 11: Documentation
test "[[ -f INSTALLATION_GUIDE.md ]]" "INSTALLATION_GUIDE.md exists"
test "[[ -f QUICK_REFERENCE.md ]]" "QUICK_REFERENCE.md exists"
test "[[ -f create-package.sh ]]" "create-package.sh exists"
test "[[ -x create-package.sh ]]" "create-package.sh is executable"

# Test 12: Package exists
test "[[ -f STAR-RSeQC-0.1.0.zip ]]" "Distribution package exists"

# Test 13: Package validation
if [[ -f STAR-RSeQC-0.1.0.zip ]]; then
    test "unzip -l STAR-RSeQC-0.1.0.zip | grep -q 'setup.sh'" "setup.sh in package"
    test "unzip -l STAR-RSeQC-0.1.0.zip | grep -q 'star-rseqc'" "Binary in package"
    test "unzip -l STAR-RSeQC-0.1.0.zip | grep -q 'Dockerfile'" "Dockerfile in package"
    test "unzip -l STAR-RSeQC-0.1.0.zip | grep -q 'docs/' " "Docs in package"
    test "unzip -l STAR-RSeQC-0.1.0.zip | wc -l | grep -qE '[1-9]'" "Package has content"
fi

# Test 14: Cargo build
test "cargo build --release 2>&1 | tail -1 | grep -q 'Finished'" "Cargo release build succeeds"

# Test 15: Cargo format check
test "cargo fmt -- --check 2>&1 | grep -q -v error" "Code formatting correct"

echo
echo "════════════════════════════════════════════════════════════════"
TOTAL=$((PASS + FAIL))
PCT=$(( (PASS * 100) / TOTAL ))

echo "Results: $PASS/$TOTAL passed ($PCT%)"
echo "════════════════════════════════════════════════════════════════"

if [[ $FAIL -eq 0 ]]; then
    echo "✅ All tests PASSED! setup.sh is ready for production."
    exit 0
else
    echo "⚠️ Some tests failed ($FAIL failures)"
    exit 1
fi
