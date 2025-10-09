steps <- c(
  "scripts/01_normalize.R",
  "scripts/02_de_limma.R",
  "scripts/03_go_enrichment.R",
  "scripts/04_plots_and_goi.R",
  "scripts/05_qc_and_volcano.R",
  "scripts/06_go_networks.R"
)

run_step <- function(s) {
  cat("\n========== RUNNING:", s, "==========\n")
  flush.console()
  tryCatch({
    # çalışma dizinini değiştirme!
    source(s, echo = FALSE, max.deparse.length = Inf, chdir = FALSE)
    cat("✅ DONE:", s, "\n")
  }, error = function(e) {
    cat("❌ FAILED at", s, ":\n", conditionMessage(e), "\n")
    quit(status = 1)
  })
}

# run from repo root
cat("Working dir:", getwd(), "\n")
invisible(lapply(steps, run_step))
cat("\n🎉 All steps finished.\n")
