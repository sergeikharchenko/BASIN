@echo off
TITLE BASIN App Launcher
ECHO Starting BASIN...
R -e "shiny::runApp('app.R', launch.browser = TRUE)"
PAUSE