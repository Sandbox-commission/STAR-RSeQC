use chrono::Local;
use crossterm::{
    cursor, execute, queue,
    style::{self, Color},
    terminal::{self, ClearType},
};
use std::io::{self, Write};
use std::time::Duration;

// ─── Public snapshot types ───────────────────────────────────────────────────

#[derive(Clone)]
pub struct JobSlotSnapshot {
    pub sample: String,
    pub step: String,
    pub elapsed_secs: f64,
    pub pct: usize,
}

#[allow(dead_code)]
pub struct RenderSnapshot {
    pub done: usize,
    pub total: usize,
    pub completed: usize,
    pub skipped: usize,
    pub failed: usize,
    pub phase_label: String,
    pub jobs: Vec<Option<JobSlotSnapshot>>,
    pub recent_events: Vec<String>,
    pub elapsed: Duration,
    pub avg_dur: f64,
    pub overall_frac: f64,
    pub overall_phase: usize,
    pub overall_total_phases: usize,
    pub overall_done: usize,
    pub overall_total: usize,
    pub overall_elapsed: Duration,
    pub p3_frac: Option<f64>,
    pub cancelled: bool,
    pub resumed: usize,
}

// ─── Layout ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum LayoutMode {
    Compact,
    Normal,
    Wide,
}

#[allow(dead_code)]
pub struct LayoutMetrics {
    pub mode: LayoutMode,
    pub overall_bar_w: usize,
    pub sample_col_w: usize,
    pub step_col_w: usize,
    pub bar_col_w: usize,
    pub pct_col_w: usize,
    pub stats_rows: usize,
}

pub fn compute_layout(w: usize, _h: usize) -> LayoutMetrics {
    let mode = match w {
        0..=79 => LayoutMode::Compact,
        80..=119 => LayoutMode::Normal,
        _ => LayoutMode::Wide,
    };

    let overall_bar_w = w.saturating_sub(7).max(4);

    // Table column widths as fractions of terminal width
    // sample 50%, step 20%, bar 20%, pct = remaining
    let sample_col_w = (w * 50 / 100).max(10);
    let step_col_w = (w * 20 / 100).max(8);
    let bar_col_w = (w * 20 / 100).max(8);
    // 5 border chars: | col | col | col | col |
    let used = sample_col_w + step_col_w + bar_col_w + 5;
    let pct_col_w = w.saturating_sub(used).max(5);

    let stats_rows = if mode == LayoutMode::Wide { 1 } else { 2 };

    LayoutMetrics {
        mode,
        overall_bar_w,
        sample_col_w,
        step_col_w,
        bar_col_w,
        pct_col_w,
        stats_rows,
    }
}

// ─── Duration formatting (pub for main.rs post-run summary) ──────────────────

pub fn fmt_duration(d: Duration) -> String {
    let secs = d.as_secs();
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    if h > 0 {
        format!("{h:02}:{m:02}:{s:02}")
    } else {
        format!("{m:02}:{s:02}")
    }
}

pub fn fmt_secs(s: f64) -> String {
    if s.is_nan() || s.is_infinite() || s < 0.0 {
        return "??:??".to_string();
    }
    if s >= 359_999.0 {
        return "99:59:59+".to_string();
    }
    fmt_duration(Duration::from_secs_f64(s))
}

pub fn truncate_to(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        s.chars().take(n).collect()
    } else {
        s.to_string()
    }
}

// ─── Gradient bar (inline ANSI) ──────────────────────────────────────────────

/// Build a gradient progress bar as a String with embedded ANSI escape codes.
fn gradient_bar_string(
    filled: usize,
    empty: usize,
    dark: (u8, u8, u8),
    bright: (u8, u8, u8),
    blink_on: bool,
    blink_color: (u8, u8, u8),
) -> String {
    let mut buf = String::with_capacity((filled + empty) * 20);
    if filled == 0 {
        buf.push_str(&fg_rgb(128, 128, 128));
        for _ in 0..empty {
            buf.push('\u{2591}');
        }
        buf.push_str(RESET);
        return buf;
    }
    for i in 0..filled {
        let t = if filled > 1 {
            i as f64 / (filled - 1) as f64
        } else {
            1.0
        };
        let is_last = i == filled - 1;
        if is_last && blink_on && empty > 0 {
            let (r, g, b) = blink_color;
            buf.push_str(&fg_rgb(r, g, b));
        } else {
            let r = (dark.0 as f64 + t * (bright.0 as f64 - dark.0 as f64)) as u8;
            let g = (dark.1 as f64 + t * (bright.1 as f64 - dark.1 as f64)) as u8;
            let b = (dark.2 as f64 + t * (bright.2 as f64 - dark.2 as f64)) as u8;
            buf.push_str(&fg_rgb(r, g, b));
        }
        buf.push('\u{2588}');
    }
    buf.push_str(&fg_rgb(128, 128, 128));
    for _ in 0..empty {
        buf.push('\u{2591}');
    }
    buf.push_str(RESET);
    buf
}

/// Set foreground color — uses true-color or 256-color depending on terminal.
fn set_fg(stdout: &mut io::Stdout, r: u8, g: u8, b: u8) {
    if is_truecolor() {
        let _ = execute!(stdout, style::SetForegroundColor(Color::Rgb { r, g, b }));
    } else {
        let _ = execute!(stdout, style::SetForegroundColor(
            Color::AnsiValue(rgb_to_256(r, g, b))
        ));
    }
}

/// Print a gradient progress bar directly to stdout (used by post-run summary in main.rs).
pub fn print_gradient_bar(
    stdout: &mut io::Stdout,
    filled: usize,
    empty: usize,
    dark: (u8, u8, u8),
    bright: (u8, u8, u8),
    blink_on: bool,
    blink_color: (u8, u8, u8),
) {
    if filled == 0 {
        let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
        for _ in 0..empty {
            let _ = write!(stdout, "\u{2591}");
        }
        let _ = stdout.flush();
        return;
    }
    for i in 0..filled {
        let t = if filled > 1 {
            i as f64 / (filled - 1) as f64
        } else {
            1.0
        };
        let is_last = i == filled - 1;
        if is_last && blink_on && empty > 0 {
            let (r, g, b) = blink_color;
            set_fg(stdout, r, g, b);
        } else {
            let r = (dark.0 as f64 + t * (bright.0 as f64 - dark.0 as f64)) as u8;
            let g = (dark.1 as f64 + t * (bright.1 as f64 - dark.1 as f64)) as u8;
            let b = (dark.2 as f64 + t * (bright.2 as f64 - dark.2 as f64)) as u8;
            set_fg(stdout, r, g, b);
        }
        let _ = write!(stdout, "\u{2588}");
    }
    let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
    for _ in 0..empty {
        let _ = write!(stdout, "\u{2591}");
    }
    let _ = stdout.flush();
}

// ─── Per-job bar color helper (matches HTML mockup) ─────────────────────────
//
// Mockup palette:
//   amber  (<25%)  : #854f0b → #ef9f27
//   teal   (25–59%): #1d9e75 → #5dcaa5
//   blue   (≥60%)  : #185fa5 → #378add

fn job_bar_colors(frac: f64) -> ((u8, u8, u8), (u8, u8, u8)) {
    if frac < 0.25 {
        // amber
        ((0x85, 0x4f, 0x0b), (0xef, 0x9f, 0x27))
    } else if frac < 0.60 {
        // teal
        ((0x1d, 0x9e, 0x75), (0x5d, 0xca, 0xa5))
    } else {
        // blue
        ((0x18, 0x5f, 0xa5), (0x37, 0x8a, 0xdd))
    }
}

/// Return the ANSI fg color string for percentage text matching the bar tier.
fn job_pct_color(frac: f64) -> String {
    if frac < 0.25 {
        fg_rgb(0xef, 0x9f, 0x27) // amber
    } else if frac < 0.60 {
        fg_rgb(0x1d, 0x9e, 0x75) // teal
    } else {
        fg_rgb(0x37, 0x8a, 0xdd) // blue
    }
}

// ─── ANSI color helpers for frame buffer lines ───────────────────────────────

use std::sync::atomic::{AtomicU8, Ordering as AtomicOrdering};

/// 0 = not checked, 1 = truecolor, 2 = 256-color
static COLOR_MODE: AtomicU8 = AtomicU8::new(0);

fn is_truecolor() -> bool {
    let v = COLOR_MODE.load(AtomicOrdering::Relaxed);
    if v != 0 {
        return v == 1;
    }
    let tc = std::env::var("COLORTERM")
        .map(|v| v == "truecolor" || v == "24bit")
        .unwrap_or(false);
    // Also check if NOT inside screen/tmux-via-screen
    let term = std::env::var("TERM").unwrap_or_default();
    let in_screen = term.starts_with("screen") && std::env::var("TMUX").is_err();
    let result = tc && !in_screen;
    COLOR_MODE.store(if result { 1 } else { 2 }, AtomicOrdering::Relaxed);
    result
}

/// Convert RGB to nearest 256-color index.
fn rgb_to_256(r: u8, g: u8, b: u8) -> u8 {
    // Check if it's close to a greyscale value (232-255: grey ramp)
    if r.abs_diff(g) < 10 && g.abs_diff(b) < 10 {
        let avg = ((r as u16 + g as u16 + b as u16) / 3) as u8;
        if avg < 8 { return 16; }
        if avg > 248 { return 231; }
        return 232 + ((avg as u16 - 8) * 24 / 240) as u8;
    }
    // Map to 6x6x6 color cube (indices 16-231)
    let ri = ((r as u16) * 5 / 255) as u8;
    let gi = ((g as u16) * 5 / 255) as u8;
    let bi = ((b as u16) * 5 / 255) as u8;
    16 + 36 * ri + 6 * gi + bi
}

/// Foreground color — true-color when supported, 256-color fallback.
fn fg_rgb(r: u8, g: u8, b: u8) -> String {
    if is_truecolor() {
        format!("\x1b[38;2;{r};{g};{b}m")
    } else {
        format!("\x1b[38;5;{}m", rgb_to_256(r, g, b))
    }
}

/// Bold + foreground color.
fn fg_rgb_bold(r: u8, g: u8, b: u8) -> String {
    if is_truecolor() {
        format!("\x1b[1;38;2;{r};{g};{b}m")
    } else {
        format!("\x1b[1;38;5;{}m", rgb_to_256(r, g, b))
    }
}

fn fg_named(name: &str) -> &'static str {
    match name {
        "white" => "\x1b[97m",
        "cyan" => "\x1b[96m",
        "green" => "\x1b[92m",
        "yellow" => "\x1b[93m",
        "red" => "\x1b[91m",
        "magenta" => "\x1b[95m",
        "darkgrey" => "\x1b[90m",
        "darkyellow" => "\x1b[33m",
        "darkred" => "\x1b[31m",
        "bold" => "\x1b[1m",
        "reset" => "\x1b[0m",
        _ => "\x1b[0m",
    }
}

const RESET: &str = "\x1b[0m";

// ─── Theme colors (resolved at runtime) ─────────────────────────────────────

/// Lazily initialized color theme using 256-color or true-color.
struct Theme {
    section_label: String,
    hdr_blue: String,
    stats_grey: String,
    stats_val: String,
    badge_blue: String,
    badge_text: String,
    sep_dim: String,
    cnt_done: String,
    cnt_skip: String,
    cnt_fail: String,
    cnt_rem: String,
    act_done: String,
    spin_color: String,
}

impl Theme {
    fn load() -> Self {
        Self {
            section_label: fg_rgb_bold(240, 230, 140),
            hdr_blue:      fg_rgb_bold(88, 166, 255),
            stats_grey:    fg_rgb_bold(139, 148, 158),
            stats_val:     fg_rgb_bold(230, 237, 243),
            badge_blue:    fg_rgb_bold(31, 111, 235),
            badge_text:    fg_rgb_bold(88, 166, 255),
            sep_dim:       fg_rgb_bold(68, 68, 68),
            cnt_done:      fg_rgb_bold(63, 185, 80),
            cnt_skip:      fg_rgb_bold(139, 148, 158),
            cnt_fail:      fg_rgb_bold(248, 81, 73),
            cnt_rem:       fg_rgb_bold(88, 166, 255),
            act_done:      fg_rgb_bold(63, 185, 80),
            spin_color:    fg_rgb_bold(239, 159, 39),
        }
    }
}

use std::sync::OnceLock;
static THEME: OnceLock<Theme> = OnceLock::new();

fn theme() -> &'static Theme {
    THEME.get_or_init(Theme::load)
}

// Accessor functions for theme colors — return &str for use in format! strings
#[allow(non_snake_case)] fn SECTION_LABEL() -> &'static str { &theme().section_label }
#[allow(non_snake_case)] fn HDR_BLUE()      -> &'static str { &theme().hdr_blue }
#[allow(non_snake_case)] fn STATS_GREY()    -> &'static str { &theme().stats_grey }
#[allow(non_snake_case)] fn STATS_VAL()     -> &'static str { &theme().stats_val }
#[allow(non_snake_case)] fn BADGE_BLUE()    -> &'static str { &theme().badge_blue }
#[allow(non_snake_case)] fn BADGE_TEXT()    -> &'static str { &theme().badge_text }
#[allow(non_snake_case)] fn SEP_DIM()       -> &'static str { &theme().sep_dim }
#[allow(non_snake_case)] fn CNT_DONE()      -> &'static str { &theme().cnt_done }
#[allow(non_snake_case)] fn CNT_SKIP()      -> &'static str { &theme().cnt_skip }
#[allow(non_snake_case)] fn CNT_FAIL()      -> &'static str { &theme().cnt_fail }
#[allow(non_snake_case)] fn CNT_REM()       -> &'static str { &theme().cnt_rem }
#[allow(non_snake_case)] fn ACT_DONE()      -> &'static str { &theme().act_done }
#[allow(non_snake_case)] fn SPIN_COLOR()    -> &'static str { &theme().spin_color }

// ─── Pad a visual line to terminal width ─────────────────────────────────────

/// Count the visible (non-ANSI) characters in a string.
fn visible_len(s: &str) -> usize {
    let mut count = 0;
    let mut in_escape = false;
    for c in s.chars() {
        if in_escape {
            if c == 'm' {
                in_escape = false;
            }
        } else if c == '\x1b' {
            in_escape = true;
        } else {
            count += 1;
        }
    }
    count
}

/// Pad line with spaces to exactly `w` visible characters to overwrite stale content.
fn pad_to_width(s: &str, w: usize) -> String {
    let vis = visible_len(s);
    if vis >= w {
        s.to_string()
    } else {
        format!("{}{}", s, " ".repeat(w - vis))
    }
}

/// Truncate a string with ANSI escapes to `max_vis` visible characters.
fn truncate_ansi(s: &str, max_vis: usize) -> String {
    let mut out = String::with_capacity(s.len());
    let mut vis = 0;
    let mut in_escape = false;
    for c in s.chars() {
        if in_escape {
            out.push(c);
            if c == 'm' {
                in_escape = false;
            }
        } else if c == '\x1b' {
            in_escape = true;
            out.push(c);
        } else {
            if vis >= max_vis {
                break;
            }
            out.push(c);
            vis += 1;
        }
    }
    out
}

// ─── Frame clipping for resize safety ────────────────────────────────────────

/// Number of rows reserved for the header section (always kept at top).
const HEADER_ROWS: usize = 4;
/// Number of rows reserved for the footer section (always kept at bottom).
const FOOTER_ROWS: usize = 3;

/// Clip frame to fit within `h` rows.
/// Keeps the first HEADER_ROWS and last FOOTER_ROWS intact.
/// Middle rows are truncated to fit. Never emits more than `h` rows.
fn clip_frame(frame: &[String], h: usize) -> Vec<String> {
    if frame.len() <= h {
        return frame.to_vec();
    }
    let header_take = HEADER_ROWS.min(frame.len()).min(h);
    let footer_take = FOOTER_ROWS.min(frame.len().saturating_sub(header_take)).min(h.saturating_sub(header_take));
    let middle_budget = h.saturating_sub(header_take + footer_take);

    let mut out = Vec::with_capacity(h);

    // Header
    for row in frame.iter().take(header_take) {
        out.push(row.clone());
    }

    // Middle: take from the section between header and footer
    let footer_start = frame.len().saturating_sub(footer_take);
    let middle_start = header_take;
    let middle_end = footer_start;
    let middle_avail = middle_end.saturating_sub(middle_start);
    if middle_budget > 0 && middle_avail > 0 {
        let take = middle_budget.min(middle_avail);
        for row in frame[middle_start..middle_start + take].iter() {
            out.push(row.clone());
        }
    }

    // Footer
    for row in frame[footer_start..].iter().take(footer_take) {
        out.push(row.clone());
    }

    // Ensure we never exceed h
    out.truncate(h);
    out
}

// ─── Braille spinners (matches HTML mockup) ─────────────────────────────────

const BRAILLE_SPINNERS: [char; 10] = [
    '\u{280B}', // ⠋
    '\u{2819}', // ⠙
    '\u{2839}', // ⠹
    '\u{2838}', // ⠸
    '\u{283C}', // ⠼
    '\u{2834}', // ⠴
    '\u{2826}', // ⠦
    '\u{2827}', // ⠧
    '\u{2807}', // ⠇
    '\u{280F}', // ⠏
];

// ─── Main render entry point ─────────────────────────────────────────────────

pub fn render(
    stdout: &mut io::Stdout,
    snap: &RenderSnapshot,
    parallel_jobs: usize,
    blink_on: bool,
) {
    // Read terminal size for frame building
    let (build_w, build_h) = terminal::size().unwrap_or((80, 24));
    let w = build_w as usize;
    let h = build_h as usize;
    if w < 20 || h < 10 {
        return;
    }

    let lm = compute_layout(w, h);

    #[allow(non_snake_case)] let SECTION_LABEL = SECTION_LABEL();
    #[allow(non_snake_case)] let HDR_BLUE = HDR_BLUE();
    #[allow(non_snake_case)] let STATS_GREY = STATS_GREY();
    #[allow(non_snake_case)] let STATS_VAL = STATS_VAL();
    #[allow(non_snake_case)] let BADGE_BLUE = BADGE_BLUE();
    #[allow(non_snake_case)] let BADGE_TEXT = BADGE_TEXT();
    #[allow(non_snake_case)] let SEP_DIM = SEP_DIM();
    #[allow(non_snake_case)] let CNT_DONE = CNT_DONE();
    #[allow(non_snake_case)] let CNT_SKIP = CNT_SKIP();
    #[allow(non_snake_case)] let CNT_FAIL = CNT_FAIL();
    #[allow(non_snake_case)] let CNT_REM = CNT_REM();
    #[allow(non_snake_case)] let ACT_DONE = ACT_DONE();
    #[allow(non_snake_case)] let SPIN_COLOR = SPIN_COLOR();

    let mut frame: Vec<String> = Vec::with_capacity(h);

    // ── Header ──
    {
        let sep = format!("{}{}{RESET}", SEP_DIM, "\u{2500}".repeat(w));
        frame.push(pad_to_width(&sep, w));
    }
    {
        let title = "STAR-RSeQC";
        let ver = env!("CARGO_PKG_VERSION");
        let title_str = format!("{title} v{ver}");
        let title_vis = title_str.len();
        let pad_l = w.saturating_sub(title_vis) / 2;
        let line = format!(
            "{}{}{}{}{RESET}",
            " ".repeat(pad_l),
            fg_named("bold"),
            HDR_BLUE,
            title_str,
        );
        frame.push(pad_to_width(&line, w));
    }
    {
        let sub = "STAR 2-pass alignment + RSeQC quality control  |  paired-end RNA-seq";
        let sub_vis = sub.len();
        let pad_l = w.saturating_sub(sub_vis) / 2;
        let line = format!(
            "{}{}{}{RESET}",
            " ".repeat(pad_l),
            STATS_GREY,
            sub,
        );
        frame.push(pad_to_width(&line, w));
    }
    {
        let sep = format!("{}{}{RESET}", SEP_DIM, "\u{2500}".repeat(w));
        frame.push(pad_to_width(&sep, w));
    }

    // ── Overall pipeline section ──
    {
        let label = format!("{}{SECTION_LABEL}OVERALL PIPELINE{RESET}", "  ");
        let phase_name = snap.phase_label.split('\u{2014}').last().unwrap_or(&snap.phase_label).trim();
        let badge = format!(
            "{BADGE_BLUE}[{BADGE_TEXT} phase {} / {} \u{2014} {} {BADGE_BLUE}]{RESET}",
            snap.overall_phase,
            snap.overall_total_phases,
            phase_name,
        );
        let label_vis = 18; // "  OVERALL PIPELINE"
        let badge_vis = visible_len(&badge);
        let gap = w.saturating_sub(label_vis + badge_vis);
        let line = format!("{label}{}{badge}", " ".repeat(gap));
        frame.push(pad_to_width(&line, w));
    }

    // Overall bar — teal gradient (#1d9e75 → #5dcaa5)
    {
        let bar_w = lm.overall_bar_w;
        let o_filled = (bar_w as f64 * snap.overall_frac) as usize;
        let o_empty = bar_w.saturating_sub(o_filled);
        let o_pct = (snap.overall_frac * 100.0) as usize;
        let bar = gradient_bar_string(
            o_filled, o_empty,
            (0x1d, 0x9e, 0x75), (0x5d, 0xca, 0xa5),
            blink_on, (0x5d, 0xca, 0xa5),
        );
        let line = format!("  {bar} {STATS_GREY}{:>3}%{RESET}", o_pct);
        frame.push(pad_to_width(&line, w));
    }

    // Overall stats row
    {
        let o_eta_str = if snap.overall_frac > 0.0 && snap.overall_frac < 1.0 {
            let o_eta = Duration::from_secs_f64(
                snap.overall_elapsed.as_secs_f64() / snap.overall_frac * (1.0 - snap.overall_frac),
            );
            fmt_duration(o_eta)
        } else if snap.overall_frac >= 1.0 {
            "00:00".to_string()
        } else {
            "--:--".to_string()
        };
        let line = format!(
            "    {STATS_GREY}elapsed {STATS_VAL}{}   {STATS_GREY}eta {STATS_VAL}{}   {STATS_GREY}samples {STATS_VAL}{}/{}{RESET}",
            fmt_duration(snap.overall_elapsed),
            o_eta_str,
            snap.overall_done, snap.overall_total,
        );
        frame.push(pad_to_width(&line, w));
    }

    // ── Phase progress section ──
    frame.push(pad_to_width(
        &format!("  {SECTION_LABEL}PHASE PROGRESS{RESET}"),
        w,
    ));

    // Phase bar computation
    let remaining = snap.total.saturating_sub(snap.done);
    let (phase_frac, phase_pct) = if let Some(p3f) = snap.p3_frac {
        (p3f, (p3f * 100.0) as usize)
    } else if snap.total > 0 {
        let f = snap.done.min(snap.total) as f64 / snap.total as f64;
        (f, (f * 100.0) as usize)
    } else {
        (0.0, 0)
    };

    // Phase bar — blue gradient (#185fa5 → #378add)
    {
        let bar_w = lm.overall_bar_w;
        let p_filled = (bar_w as f64 * phase_frac) as usize;
        let p_empty = bar_w.saturating_sub(p_filled);
        let bar = gradient_bar_string(
            p_filled, p_empty,
            (0x18, 0x5f, 0xa5), (0x37, 0x8a, 0xdd),
            blink_on, (0x37, 0x8a, 0xdd),
        );
        let line = format!("  {bar} {STATS_GREY}{:>3}%{RESET}", phase_pct);
        frame.push(pad_to_width(&line, w));
    }

    // Phase stats
    let avg_dur = snap.avg_dur;
    let processed = snap.completed;
    let elapsed = snap.elapsed;
    let elapsed_mins = elapsed.as_secs_f64() / 60.0;
    let speed_per_min = if elapsed_mins > 0.01 && processed > 0 {
        processed as f64 / elapsed_mins
    } else {
        0.0
    };
    // Show /hr for slow jobs (< 1/min), "calculating..." before first completion
    let speed_str = if processed == 0 {
        "calculating...".to_string()
    } else if speed_per_min < 1.0 {
        format!("{:.1}/hr", speed_per_min * 60.0)
    } else {
        format!("{:.1}/min", speed_per_min)
    };
    let has_eta;
    let eta = if processed > 0 && remaining > 0 {
        has_eta = true;
        Duration::from_secs_f64(avg_dur * remaining as f64)
    } else if snap.done > 0 && remaining > 0 {
        has_eta = true;
        let elapsed_secs = elapsed.as_secs_f64().max(0.001);
        let per = elapsed_secs / snap.done as f64;
        Duration::from_secs_f64(per * remaining as f64)
    } else {
        has_eta = false;
        Duration::ZERO
    };
    let eta_str = if has_eta {
        fmt_duration(eta)
    } else {
        "--:--".to_string()
    };

    let completion_time = if has_eta && eta.as_secs() > 0 {
        let now = Local::now();
        let duration_secs = eta.as_secs() as i64;
        let completion = now + chrono::Duration::seconds(duration_secs);
        completion.format("%H:%M:%S").to_string()
    } else {
        "\u{2014}".to_string()
    };

    if lm.stats_rows == 2 {
        let s1 = format!(
            "    {STATS_GREY}elapsed {STATS_VAL}{}   {STATS_GREY}eta {STATS_VAL}{}   {STATS_GREY}{}/{} done{RESET}",
            fmt_duration(elapsed),
            eta_str,
            snap.done, snap.total,
        );
        frame.push(pad_to_width(&s1, w));
        let s2 = format!(
            "    {STATS_GREY}speed {STATS_VAL}{}   {STATS_GREY}complete {STATS_VAL}{}{RESET}",
            speed_str, completion_time,
        );
        frame.push(pad_to_width(&s2, w));
    } else {
        let s = format!(
            "    {STATS_GREY}elapsed {STATS_VAL}{}   {STATS_GREY}eta {STATS_VAL}{}   {STATS_GREY}speed {STATS_VAL}{}   {STATS_GREY}complete {STATS_VAL}{}{RESET}",
            fmt_duration(elapsed),
            eta_str, speed_str, completion_time,
        );
        frame.push(pad_to_width(&s, w));
    }

    // ── Active jobs section ──
    let active_count = snap.jobs.iter().filter(|s| s.is_some()).count();
    frame.push(pad_to_width(
        &format!(
            "  {SECTION_LABEL}ACTIVE JOBS {}{CNT_REM}({}/{}){RESET}",
            RESET, active_count, parallel_jobs,
        ),
        w,
    ));

    let spin_idx = (snap.elapsed.as_millis() / 100) as usize;

    let sc = lm.sample_col_w;
    let stc = lm.step_col_w;
    let bc = lm.bar_col_w;
    let pc = lm.pct_col_w;

    // Table header
    {
        let hdr = format!(
            "{SEP_DIM}\u{250C}{}\u{252C}{}\u{252C}{}\u{252C}{}\u{2510}{RESET}",
            "\u{2500}".repeat(sc),
            "\u{2500}".repeat(stc),
            "\u{2500}".repeat(bc),
            "\u{2500}".repeat(pc),
        );
        frame.push(pad_to_width(&hdr, w));
    }
    {
        let sample_hdr = format!("{:<width$}", " SAMPLE", width = sc);
        let step_hdr = format!("{:<width$}", " STEP", width = stc);
        let prog_hdr = format!("{:<width$}", " PROGRESS", width = bc);
        let pct_hdr = format!("{:>width$}", " ", width = pc);
        let line = format!(
            "{SEP_DIM}\u{2502}{STATS_GREY}{sample_hdr}{SEP_DIM}\u{2502}{STATS_GREY}{step_hdr}{SEP_DIM}\u{2502}{STATS_GREY}{prog_hdr}{SEP_DIM}\u{2502}{STATS_GREY}{pct_hdr}{SEP_DIM}\u{2502}{RESET}",
        );
        frame.push(pad_to_width(&line, w));
    }
    {
        let sep = format!(
            "{SEP_DIM}\u{251C}{}\u{253C}{}\u{253C}{}\u{253C}{}\u{2524}{RESET}",
            "\u{2500}".repeat(sc),
            "\u{2500}".repeat(stc),
            "\u{2500}".repeat(bc),
            "\u{2500}".repeat(pc),
        );
        frame.push(pad_to_width(&sep, w));
    }

    // Job rows
    let active_jobs: Vec<(usize, &JobSlotSnapshot)> = snap
        .jobs
        .iter()
        .enumerate()
        .filter_map(|(i, s)| s.as_ref().map(|j| (i, j)))
        .collect();

    if active_jobs.is_empty() {
        let empty_cell = format!("{:<width$}", " No active jobs", width = sc + stc + bc + pc + 3);
        let line = format!(
            "{SEP_DIM}\u{2502}{STATS_GREY}{empty_cell}{SEP_DIM}\u{2502}{RESET}",
        );
        frame.push(pad_to_width(&line, w));
    }

    // Limit jobs shown based on terminal height.
    let rows_below_jobs = 7;
    let max_job_pairs = h.saturating_sub(frame.len() + rows_below_jobs) / 2;
    let max_job_rows = max_job_pairs.max(1);

    for (shown, (i, job)) in active_jobs.iter().enumerate() {
        if shown >= max_job_rows {
            let hidden = active_jobs.len().saturating_sub(shown);
            if hidden > 0 {
                let more_cell = format!(
                    "{:<width$}",
                    format!(" ... and {hidden} more"),
                    width = sc + stc + bc + pc + 3
                );
                let line = format!(
                    "{SEP_DIM}\u{2502}{STATS_GREY}{more_cell}{SEP_DIM}\u{2502}{RESET}",
                );
                frame.push(pad_to_width(&line, w));
            }
            break;
        }

        let spin = BRAILLE_SPINNERS[(spin_idx + i) % BRAILLE_SPINNERS.len()];

        // Sample cell
        let name_max = sc.saturating_sub(4); // " S name "
        let name_display = if job.sample.chars().count() > name_max {
            format!(
                "{}...",
                job.sample
                    .chars()
                    .take(name_max.saturating_sub(3))
                    .collect::<String>()
            )
        } else {
            job.sample.clone()
        };
        let sample_cell = format!("{:<width$}", format!(" {spin} {name_display}"), width = sc);

        // Step cell
        let step_display = truncate_to(&job.step, stc.saturating_sub(2));
        let step_cell = format!("{:<width$}", format!(" {step_display}"), width = stc);

        // Determine fill fraction
        let (frac, pct_display, is_real_pct) = if job.pct > 0 {
            let f = job.pct.min(100) as f64 / 100.0;
            (f, format!("{:>3}%", job.pct.min(100)), true)
        } else if avg_dur > 0.0 {
            let f = (job.elapsed_secs / avg_dur).min(1.0);
            (f, format!("~{:>2}%", (f * 100.0) as usize), false)
        } else {
            (0.0, String::new(), false)
        };

        // Bar cell
        let bar_inner_w = bc.saturating_sub(4); // " [bar] "
        let bar_cell = if avg_dur > 0.0 || is_real_pct {
            let b_filled = (bar_inner_w as f64 * frac) as usize;
            let b_empty = bar_inner_w.saturating_sub(b_filled);
            let (dark, bright) = job_bar_colors(frac);
            let bar = gradient_bar_string(b_filled, b_empty, dark, bright, blink_on, bright);
            format!(" [{}] ", bar)
        } else {
            // Pulse animation when no avg_dur
            let pulse_pos = (spin_idx + i * 3) % (bar_inner_w + 4);
            let mut pulse = String::new();
            for p in 0..bar_inner_w {
                if p >= pulse_pos.saturating_sub(2) && p <= pulse_pos {
                    pulse.push_str(&format!("{}\u{2588}", fg_rgb(0xef, 0x9f, 0x27)));
                } else {
                    pulse.push_str(&format!("{}\u{2591}", fg_rgb(128, 128, 128)));
                }
            }
            format!(" [{}] ", pulse)
        };
        // Pad bar_cell visually to bc
        let bar_cell_vis = visible_len(&bar_cell);
        let bar_cell_padded = if bar_cell_vis < bc {
            format!("{}{}", bar_cell, " ".repeat(bc - bar_cell_vis))
        } else {
            bar_cell
        };

        // Pct cell — color matches bar tier
        let pct_color = job_pct_color(frac);
        let pct_cell = format!("{pct_color}{:>width$}{RESET}", pct_display, width = pc);

        let line = format!(
            "{SEP_DIM}\u{2502}{SPIN_COLOR}{sample_cell}{SEP_DIM}\u{2502}{}{step_cell}{SEP_DIM}\u{2502}{}{bar_cell_padded}{SEP_DIM}\u{2502}{pct_cell}{SEP_DIM}\u{2502}{RESET}",
            HDR_BLUE, "",
        );
        frame.push(pad_to_width(&line, w));

        // Second row: elapsed / ~ETA (spans full table width)
        let ela_str = fmt_secs(job.elapsed_secs);
        let eta_part = if job.pct > 0 && job.pct < 100 {
            let total_est = job.elapsed_secs / (job.pct as f64 / 100.0);
            let rem = total_est - job.elapsed_secs;
            format!("{} / ~{}", ela_str, fmt_secs(rem.max(0.0)))
        } else if avg_dur > 0.0 {
            format!("{} / ~{}", ela_str, fmt_secs(avg_dur))
        } else {
            ela_str
        };
        let inner_w = sc + stc + bc + pc + 3;
        let elapsed_cell = format!("{:<width$}", format!("    {eta_part}"), width = inner_w);
        let line2 = format!(
            "{SEP_DIM}\u{2502}{STATS_GREY}{elapsed_cell}{SEP_DIM}\u{2502}{RESET}",
        );
        frame.push(pad_to_width(&line2, w));
    }

    // Table bottom
    {
        let bot = format!(
            "{SEP_DIM}\u{2514}{}\u{2534}{}\u{2534}{}\u{2534}{}\u{2518}{RESET}",
            "\u{2500}".repeat(sc),
            "\u{2500}".repeat(stc),
            "\u{2500}".repeat(bc),
            "\u{2500}".repeat(pc),
        );
        frame.push(pad_to_width(&bot, w));
    }

    // ── Counters (matches HTML .summary-bar) ──
    {
        let sep = format!("{}{}{RESET}", SEP_DIM, "\u{2500}".repeat(w));
        frame.push(pad_to_width(&sep, w));
    }
    {
        let line = format!(
            "  {CNT_DONE}\u{2713} completed: {}   {CNT_SKIP}\u{2192} skipped: {}   {CNT_FAIL}\u{2717} failed: {}   {CNT_REM}\u{22C5} remaining: {}{RESET}",
            snap.completed,
            snap.skipped,
            snap.failed,
            remaining,
        );
        frame.push(pad_to_width(&line, w));
    }

    // ── Recent activity ──
    {
        let sep = format!("{}{}{RESET}", SEP_DIM, "\u{2500}".repeat(w));
        frame.push(pad_to_width(&sep, w));
    }
    frame.push(pad_to_width(
        &format!("  {SECTION_LABEL}RECENT ACTIVITY{RESET}"),
        w,
    ));

    // Leave room for footer (2 rows: separator + footer)
    let used_rows = frame.len();
    let max_event_rows = h.saturating_sub(used_rows + 2);
    let events = &snap.recent_events;
    let start = events.len().saturating_sub(max_event_rows);
    for ev_line in events[start..].iter().rev() {
        let ev = truncate_to(ev_line, w);
        // Parse event line to apply mockup colors: "  DONE  sample — task"
        let color = if ev.contains("DONE") {
            ACT_DONE
        } else if ev.contains("SKIP") || ev.contains("RESUME") {
            CNT_SKIP
        } else if ev.contains("FAIL") {
            CNT_FAIL
        } else if ev.contains("STOP") {
            fg_named("darkred")
        } else if ev.contains("INFO") {
            STATS_GREY
        } else {
            STATS_VAL
        };
        frame.push(pad_to_width(&format!("{color}{ev}{RESET}"), w));
    }

    // ── Fill remaining rows up to footer ──
    while frame.len() < h.saturating_sub(2) {
        frame.push(" ".repeat(w));
    }

    // ── Footer separator ──
    frame.push(pad_to_width(
        &format!("{}{}{RESET}", SEP_DIM, "\u{2500}".repeat(w)),
        w,
    ));

    // ── Footer ──
    {
        let quit_hint = if snap.cancelled {
            format!("{CNT_FAIL}  CANCELLING...")
        } else {
            format!("{STATS_GREY}  [q] quit   Ctrl+C cancel")
        };
        let timestamp = format!("Updated: {}", Local::now().format("%H:%M:%S"));
        let quit_vis = if snap.cancelled { 16 } else { 28 };
        let ts_len = timestamp.len();
        let pad = w.saturating_sub(quit_vis + ts_len);
        let line = format!(
            "{}{}{STATS_GREY}{}{RESET}",
            quit_hint,
            " ".repeat(pad),
            timestamp,
        );
        frame.push(pad_to_width(&line, w));
    }

    // ── Flush with resize safety ──
    let (flush_w, flush_h) = terminal::size().unwrap_or((build_w, build_h));
    let fw = flush_w as usize;
    let fh = flush_h as usize;

    // Clip frame to fit within fresh terminal height
    let clipped = if frame.len() > fh {
        clip_frame(&frame, fh)
    } else {
        frame
    };

    // Build entire frame into a single buffer, then one write + flush.
    // Each line is positioned absolutely with MoveTo to avoid scroll artifacts.
    let mut buf = Vec::with_capacity(fw * fh * 4);
    let _ = queue!(buf, cursor::Hide);
    for (row, line) in clipped.iter().enumerate() {
        let _ = queue!(buf, cursor::MoveTo(0, row as u16));
        let truncated = truncate_ansi(line, fw);
        let padded = pad_to_width(&truncated, fw);
        let _ = write!(buf, "{padded}");
    }
    // Clear any leftover lines below the frame
    for row in clipped.len()..fh {
        let _ = queue!(buf, cursor::MoveTo(0, row as u16));
        let _ = queue!(buf, terminal::Clear(ClearType::CurrentLine));
    }
    let _ = queue!(buf, style::ResetColor, cursor::Hide);
    // Single atomic write to stdout
    let _ = stdout.write_all(&buf);
    let _ = stdout.flush();
}

// ─── Post-run final frame ────────────────────────────────────────────────────

pub fn render_final_frame(stdout: &mut io::Stdout, message: &str) {
    let (tw, th) = terminal::size().unwrap_or((80, 24));
    let w = tw as usize;
    let h = th as usize;

    #[allow(non_snake_case)] let CNT_DONE = CNT_DONE();
    #[allow(non_snake_case)] let STATS_GREY = STATS_GREY();
    #[allow(non_snake_case)] let STATS_VAL = STATS_VAL();

    let mut frame: Vec<String> = Vec::with_capacity(h);

    // Top border
    frame.push(pad_to_width(
        &format!("{}{}{RESET}", CNT_DONE, "\u{2550}".repeat(w)),
        w,
    ));

    // Blank lines to center
    let mid = h / 2;
    for _ in 1..mid.saturating_sub(1) {
        frame.push(" ".repeat(w));
    }

    // DONE banner
    let banner_text = format!("  \u{2714} PIPELINE COMPLETE  ");
    let banner_vis = banner_text.chars().count();
    let pad_l = w.saturating_sub(banner_vis) / 2;
    let banner_line = format!(
        "{}{}{}{}{RESET}",
        " ".repeat(pad_l),
        fg_named("bold"),
        CNT_DONE,
        banner_text,
    );
    frame.push(pad_to_width(&banner_line, w));

    // Message
    let msg_vis = message.chars().count();
    let msg_pad = w.saturating_sub(msg_vis) / 2;
    let msg_line = format!(
        "{}{}{}{RESET}",
        " ".repeat(msg_pad),
        STATS_VAL,
        message,
    );
    frame.push(pad_to_width(&msg_line, w));

    // Hint
    let hint = "Press any key or wait 10s...";
    let hint_pad = w.saturating_sub(hint.len()) / 2;
    let hint_line = format!(
        "{}{}{}{RESET}",
        " ".repeat(hint_pad),
        STATS_GREY,
        hint,
    );
    frame.push(pad_to_width(&hint_line, w));

    // Fill rest
    while frame.len() < h.saturating_sub(1) {
        frame.push(" ".repeat(w));
    }

    // Bottom border
    frame.push(pad_to_width(
        &format!("{}{}{RESET}", CNT_DONE, "\u{2550}".repeat(w)),
        w,
    ));

    // Single buffered flush
    let mut buf = Vec::with_capacity(w * h * 4);
    let _ = queue!(buf, cursor::MoveTo(0, 0));
    for line in frame.iter().take(h) {
        let truncated = truncate_ansi(line, w);
        let padded = pad_to_width(&truncated, w);
        let _ = write!(buf, "{padded}\r\n");
    }
    let _ = queue!(buf, terminal::Clear(ClearType::FromCursorDown));
    let _ = queue!(buf, style::ResetColor);
    let _ = stdout.write_all(&buf);
    let _ = stdout.flush();
}
