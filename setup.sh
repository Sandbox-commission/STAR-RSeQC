#!/bin/bash
set -euo pipefail

# STAR-RSeQC Setup Script — Enhanced Version
# Three installation paths:
# 1. Docker containerization (automated)
# 2. Automatic mamba/conda install + environment setup
# 3. Manual installation with web links

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ─── Helper Functions ────────────────────────────────────────────────────────

log_header() {
    echo -e "\n${BLUE}=== $1 ===${NC}\n"
}

log_step() {
    echo -e "${YELLOW}[${1}]${NC} $2"
}

log_success() {
    echo -e "${GREEN}✓${NC} $1"
}

log_error() {
    echo -e "${RED}✗${NC} $1"
}

pause_and_ask() {
    read -p "$(echo -e ${YELLOW})Press Enter to continue...$(echo -e ${NC})" -r
}

# ─── Main Menu ──────────────────────────────────────────────────────────────

log_header "STAR-RSeQC Setup"
echo "Choose your installation method:"
echo
echo "  1) Docker Container (Recommended) — everything pre-configured"
echo "  2) Mamba/Conda (Automatic) — auto-install dependencies"
echo "  3) Manual Installation — follow web links & instructions"
echo "  4) Exit"
echo
read -p "Select option [1-4]: " -r CHOICE

case "$CHOICE" in
    1) setup_docker ;;
    2) setup_mamba_auto ;;
    3) setup_manual ;;
    4) echo "Exiting."; exit 0 ;;
    *) log_error "Invalid choice."; exit 1 ;;
esac

# ─── SECTION 1: DOCKER SETUP ────────────────────────────────────────────────

setup_docker() {
    log_header "Docker Container Installation"

    # Check if Docker is installed
    if command -v docker &> /dev/null; then
        log_success "Docker found: $(docker --version)"
        docker_setup_existing
    else
        log_error "Docker is not installed."
        echo
        echo "  a) Auto-install Docker"
        echo "  b) Manual installation (with links)"
        echo "  c) Go back to main menu"
        echo
        read -p "Select option [a-c]: " -r DOCKER_CHOICE
        case "$DOCKER_CHOICE" in
            a) install_docker_auto ;;
            b) install_docker_manual ;;
            c) main ;;
            *) log_error "Invalid choice."; exit 1 ;;
        esac
    fi
}

install_docker_auto() {
    log_step "1/4" "Detecting Linux distribution..."

    if [[ -f /etc/os-release ]]; then
        . /etc/os-release
        OS=$ID
    else
        log_error "Could not detect OS"
        exit 1
    fi

    log_step "2/4" "Installing Docker prerequisites..."
    case "$OS" in
        ubuntu|debian)
            sudo apt-get update
            sudo apt-get install -y \
                apt-transport-https \
                ca-certificates \
                curl \
                gnupg \
                lsb-release
            ;;
        fedora|rhel|centos)
            sudo dnf install -y \
                curl \
                gnupg \
                lsb-release
            ;;
        *)
            log_error "Unsupported OS: $OS"
            echo "Please visit: https://docs.docker.com/engine/install/"
            exit 1
            ;;
    esac

    log_step "3/4" "Installing Docker..."
    case "$OS" in
        ubuntu|debian)
            curl -fsSL https://download.docker.com/linux/${OS}/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
            echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/${OS} $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
            sudo apt-get update
            sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin
            ;;
        fedora|rhel|centos)
            sudo dnf config-manager --add-repo https://download.docker.com/linux/fedora/docker-ce.repo
            sudo dnf install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin
            ;;
    esac

    log_step "4/4" "Starting Docker service..."
    sudo systemctl start docker
    sudo systemctl enable docker
    sudo usermod -aG docker "$USER" || true

    log_success "Docker installed!"
    echo
    echo -e "${YELLOW}Note:${NC} You may need to run the following to use Docker without sudo:"
    echo "  newgrp docker"
    echo
    docker_setup_existing
}

install_docker_manual() {
    log_header "Manual Docker Installation"
    echo
    echo "Visit the official Docker installation guide:"
    echo "  ${BLUE}https://docs.docker.com/engine/install/${NC}"
    echo
    echo "Installation steps:"
    echo "  1. Choose your operating system"
    echo "  2. Follow the installation instructions"
    echo "  3. Verify: docker --version"
    echo "  4. Re-run this setup script after installation"
    echo
    pause_and_ask
}

docker_setup_existing() {
    log_step "1/4" "Checking docker-compose..."
    if ! command -v docker-compose &> /dev/null && ! docker compose version &> /dev/null; then
        log_error "docker-compose not found"
        exit 1
    fi
    log_success "docker-compose is available"

    log_step "2/4" "Building Docker image..."
    docker build -t star-rseqc:latest .
    log_success "Docker image built"

    log_step "3/4" "Creating config directory..."
    mkdir -p ~/.config/star-rseqc
    log_success "Config directory created"

    log_step "4/4" "Generating docker-compose command..."
    cat > docker-run.sh << 'DOCKER_SCRIPT'
#!/bin/bash
# STAR-RSeQC Docker Runner

FASTQ_DIR="${1:-.}"
OUTPUT_DIR="${2:-./star-rseqc-results}"
GENOME_DIR="${GENOME_DIR:-.}"
GTF_DIR="${GTF_DIR:-.}"

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "Error: FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Running STAR-RSeQC in Docker..."
echo "  Input:  $FASTQ_DIR"
echo "  Output: $OUTPUT_DIR"
echo

docker run --rm \
    -v "$FASTQ_DIR:/data/input:ro" \
    -v "$OUTPUT_DIR:/data/output" \
    -v "$HOME/.config/star-rseqc:/root/.config/star-rseqc:ro" \
    --cpus=16 \
    --memory=32g \
    star-rseqc:latest \
    star-rseqc /data/input -o /data/output "$@"
DOCKER_SCRIPT
    chmod +x docker-run.sh
    log_success "Created docker-run.sh"

    echo
    echo -e "${GREEN}✓ Docker setup complete!${NC}"
    echo
    echo "Usage:"
    echo "  ${YELLOW}./docker-run.sh /path/to/fastq [output_dir]${NC}"
    echo
    echo "Or with docker-compose:"
    echo "  ${YELLOW}FASTQ_DIR=/path/to/fastq docker-compose up${NC}"
    echo
}

# ─── SECTION 2: MAMBA AUTOMATIC SETUP ────────────────────────────────────────

setup_mamba_auto() {
    log_header "Automatic Mamba Installation"

    # Check if mamba/conda exists
    CONDA_ROOT=""
    for dir in ~/miniforge3 ~/mambaforge ~/miniconda3; do
        if [[ -d "$dir" ]]; then
            CONDA_ROOT="$dir"
            log_success "Found conda at: $CONDA_ROOT"
            setup_environments "$CONDA_ROOT"
            return
        fi
    done

    if [[ -d /opt/conda ]]; then
        CONDA_ROOT="/opt/conda"
        log_success "Found conda at: $CONDA_ROOT"
        setup_environments "$CONDA_ROOT"
        return
    fi

    # Conda not found — ask user
    log_error "Conda/Mamba not found"
    echo
    read -p "Auto-install Miniforge (mamba)? [Y/n] " -r INSTALL_MAMBA
    if [[ $INSTALL_MAMBA =~ ^[Yy]?$ ]]; then
        install_miniforge_auto
    else
        log_error "Mamba required for automatic setup. Please install manually:"
        install_conda_manual
        exit 1
    fi
}

install_miniforge_auto() {
    log_step "1/3" "Downloading Miniforge..."
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    MINIFORGE_INSTALLER="/tmp/Miniforge3-Linux-x86_64.sh"

    if ! curl -fsSL "$MINIFORGE_URL" -o "$MINIFORGE_INSTALLER"; then
        log_error "Failed to download Miniforge"
        exit 1
    fi
    log_success "Downloaded"

    log_step "2/3" "Installing Miniforge..."
    CONDA_ROOT="$HOME/miniforge3"
    bash "$MINIFORGE_INSTALLER" -b -p "$CONDA_ROOT"
    rm -f "$MINIFORGE_INSTALLER"
    log_success "Installed to $CONDA_ROOT"

    log_step "3/3" "Initializing conda..."
    "$CONDA_ROOT/bin/conda" init bash
    log_success "Conda initialized"

    echo
    echo -e "${YELLOW}Note:${NC} Run ${BLUE}source ~/.bashrc${NC} to activate conda"
    echo
    setup_environments "$CONDA_ROOT"
}

install_conda_manual() {
    log_header "Manual Conda Installation"
    echo
    echo "Choose your preferred conda distribution:"
    echo
    echo "  ${BLUE}Miniforge (Recommended)${NC}"
    echo "    https://github.com/conda-forge/miniforge"
    echo
    echo "  ${BLUE}Mambaforge${NC}"
    echo "    https://github.com/conda-forge/miniforge"
    echo
    echo "  ${BLUE}Miniconda${NC}"
    echo "    https://docs.conda.io/en/latest/miniconda.html"
    echo
    echo "Installation steps:"
    echo "  1. Download the appropriate installer for your system"
    echo "  2. Run: bash ~/Downloads/Miniforge3-Linux-x86_64.sh"
    echo "  3. Follow the prompts"
    echo "  4. Run: source ~/.bashrc"
    echo "  5. Re-run this setup script"
    echo
    pause_and_ask
}

setup_environments() {
    local CONDA_ROOT="$1"
    log_step "1/6" "Creating STAR environment..."
    "$CONDA_ROOT/bin/mamba" create -n star -c bioconda -c conda-forge star=2.7.11b samtools -y || true
    log_success "STAR environment ready"

    log_step "2/6" "Creating RSeQC environment..."
    "$CONDA_ROOT/bin/mamba" create -n rseqc -c bioconda -c conda-forge rseqc python -y || true
    log_success "RSeQC environment ready"

    log_step "3/6" "Finding samtools..."
    SAMTOOLS="$CONDA_ROOT/envs/star/bin/samtools"
    if [[ ! -x "$SAMTOOLS" ]]; then
        SAMTOOLS=$(command -v samtools || echo "")
    fi
    log_success "samtools: $SAMTOOLS"

    log_step "4/6" "Prompting for reference files..."
    read -p "Path to STAR genome index: " -r GENOME_DIR
    GENOME_DIR="${GENOME_DIR/#\~/$HOME}"

    read -p "Path to GTF annotation: " -r GTF_FILE
    GTF_FILE="${GTF_FILE/#\~/$HOME}"

    log_step "5/6" "Building binary..."
    cargo build --release
    mkdir -p ~/.local/bin
    cp target/release/star-rseqc ~/.local/bin/
    log_success "Binary installed to ~/.local/bin/star-rseqc"

    log_step "6/6" "Writing configuration..."
    mkdir -p ~/.config/star-rseqc
    cat > ~/.config/star-rseqc/config.json << EOF
{
  "genome_dir": "$GENOME_DIR",
  "gtf": "$GTF_FILE",
  "star_env": "$CONDA_ROOT/envs/star",
  "rseqc_env": "$CONDA_ROOT/envs/rseqc",
  "samtools": "$SAMTOOLS"
}
EOF
    log_success "Configuration written"

    echo
    echo -e "${GREEN}✓ Mamba setup complete!${NC}"
    echo
    echo "Add to PATH:"
    echo "  ${YELLOW}export PATH=\"\$HOME/.local/bin:\$PATH\"${NC}"
    echo
    echo "Run STAR-RSeQC:"
    echo "  ${YELLOW}star-rseqc /path/to/fastq${NC}"
    echo
}

# ─── SECTION 3: MANUAL SETUP ────────────────────────────────────────────────

setup_manual() {
    log_header "Manual Installation Guide"
    echo
    echo "Follow these steps to set up STAR-RSeQC manually:"
    echo
    echo -e "${BLUE}Step 1: Install Conda/Mamba${NC}"
    echo "  https://github.com/conda-forge/miniforge"
    echo "  bash ~/Downloads/Miniforge3-Linux-x86_64.sh"
    echo
    echo -e "${BLUE}Step 2: Create Environments${NC}"
    echo "  mamba create -n star -c bioconda -c conda-forge star=2.7.11b samtools"
    echo "  mamba create -n rseqc -c bioconda -c conda-forge rseqc python"
    echo
    echo -e "${BLUE}Step 3: Build STAR-RSeQC${NC}"
    echo "  git clone https://github.com/Sandbox-commission/STAR-RSeQC.git"
    echo "  cd STAR-RSeQC"
    echo "  cargo build --release"
    echo "  cp target/release/star-rseqc ~/.local/bin/"
    echo
    echo -e "${BLUE}Step 4: Create Config${NC}"
    echo "  mkdir -p ~/.config/star-rseqc"
    echo "  Create ~/.config/star-rseqc/config.json with your paths"
    echo
    echo -e "${BLUE}Step 5: Add to PATH${NC}"
    echo "  export PATH=\"\$HOME/.local/bin:\$PATH\""
    echo
    echo -e "${BLUE}Step 6: Run${NC}"
    echo "  star-rseqc /path/to/fastq"
    echo
    pause_and_ask
}

# ─── Execute based on shell mode ─────────────────────────────────────────────

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Script is being executed directly, not sourced
    true
else
    # Script is being sourced
    true
fi
