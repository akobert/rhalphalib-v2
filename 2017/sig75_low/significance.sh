echo "Plotting Significance"
combine -M Significance --signif output/testModel/model_combined.root --cminDefaultMinimizerStrategy 0 -t -1 -v 9 --expectSignal 1
