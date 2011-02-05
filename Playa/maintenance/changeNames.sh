#!/bin/bash


for f in $@; do

echo $f
#sed -i 's/SundanceExceptions.hpp/PlayaExceptions.hpp/g' $f
#sed -i 's/SundanceTabs.hpp/PlayaTabs.hpp/g' $f
#sed -i 's/SundanceHandleable.hpp/PlayaHandleable.hpp/g' $f
#sed -i 's/SundanceHandle.hpp/PlayaHandle.hpp/g' $f
#sed -i 's/SundancePrintable.hpp/PlayaPrintable.hpp/g' $f
#sed -i 's/TSFExtended/Playa/g' $f
#sed -i 's/setVerbosity/setVerb/g' $f
#sed -i 's/TSF/Playa/g' $f
#sed -i 's/ProductVector/BlockVector/g' $f
#sed -i 's/productSpace/blockSpace/g' $f
#sed -i 's/RuntimeError/std::runtime_error/g' $f
#sed -i 's/InternalError/std::logic_error/g' $f
#sed -i 's/Playa::Tabs t/Tabs t/g' $f
#sed -i 's/Sundance::Handle/Playa::Handle/g' $f
#sed -i 's/Sundance::Printable/Playa::Printable/g' $f
#sed -i 's/public Printable/public Playa::Printable/g' $f
#sed -i 's/const Printable/const Playa::Printable/g' $f

sed -i 's/NOX_Playa_Group.H/NOX_Playa_Group.hpp/g' $f
sed -i 's/PlayaNOXSolver.H/PlayaNOXSolver.hpp/g' $f

done

