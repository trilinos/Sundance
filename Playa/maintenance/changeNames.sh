#!/bin/sh

sed -i 's/SundanceOut.hpp/PlayaOut.hpp/' $0
sed -i 's/SundanceExceptions.hpp/PlayaExceptions.hpp/' $0
sed -i 's/SundanceTabs.hpp/PlayaTabs.hpp/' $0
sed -i 's/SundanceHandleable.hpp/PlayaHandleable.hpp/' $0
sed -i 's/SundanceHandle.hpp/PlayaHandle.hpp/' $0
sed -i 's/SundancePrintable.hpp/PlayaPrintable.hpp/' $0
sed -i 's/TSFExtended/Playa' $0
