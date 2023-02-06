#!/bin/bash

parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Academ_YW/ ::: data/Repbase/Academ_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/CMC_YW/ ::: data/Repbase/CMC_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Ginger_YW/ ::: data/Repbase/Ginger_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/hAT_YW/ ::: data/Repbase/hAT_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Kobolok_YW/ ::: data/Repbase/Kobolok_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Merlin_YW/ ::: data/Repbase/Merlin_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/MuDR_YW/ ::: data/Repbase/MuDR_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Novosib_YW/ ::: data/Repbase/Novosib_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/PHIS_YW/ ::: data/Repbase/PHIS_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/piggyBac_YW/ ::: data/Repbase/piggyBac_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/P_YW/ ::: data/Repbase/P_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Tc1marPlm_YW/ ::: data/Repbase/Tc1marPlm_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Transib_YW/ ::: data/Repbase/Transib_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Zator_YW/ ::: data/Repbase/Zator_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Sola1_YW/ ::: data/Repbase/Sola1_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Sola2_YW/ ::: data/Repbase/Sola2_YW/xml/*xml
parallel -j 128 --bar ./scripts/repbase_xml_parser.py -i {} -o data/Repbase/Sola3_YW/ ::: data/Repbase/Sola3_YW/xml/*xml
