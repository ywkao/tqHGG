#!/bin/bash

echo "grep -ic error plots_hadronic/log/info_stack_plots_hadronic:" `grep -ic error plots_hadronic/log/info_stack_plots_hadronic`
echo "grep -ic error plots_leptonic/log/info_stack_plots_leptonic:" `grep -ic error plots_leptonic/log/info_stack_plots_leptonic`
grep -ic error log/std*
echo ""
