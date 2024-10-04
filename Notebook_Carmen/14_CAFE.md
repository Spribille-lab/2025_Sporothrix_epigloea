CAFE5 analysis

* Coerced phylogenomic tree into ultrametric tree

using `scripts/14_CAFE.R`

* processed dbCAN4 results into format readable by CAFE5

using  `scripts/10_process_dbcanV12.R`

* first test for the number of discrete rate categories (*K*) produces the greatest likelihood

```
category=cazymefamily

cd /data/ccallen/2023_02_Sporothrix/14_CAFE
tree=/data/ccallen/2023_02_Sporothrix/09_orthology/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree
input=/data/ccallen/2023_02_Sporothrix/10_annotations/dbcan4/summarized_outputs/cazy_summarized_CAFE.txt

for n in {1..10}
do
output="$category"_singlelambda_k"$n"
mkdir "$output"
/data/ccallen/bin/CAFE5/bin/cafe5 -t "$tree" -i "$input" -p -k "$n" -o "$output" &>> "$output"/cafe.log
done

# compile results of search for optimum k
for n in {1..10}
do
echo "$category"_singlelambda_k"$n" >> "$category"_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"/Base_results.txt >> "$category"_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"/Gamma_results.txt >> "$category"_discrete_rate_categories.txt
done

# count significant families at the p=0.05 threshold:
for n in {1..10}
do
echo "$category"_singlelambda_k"$n" >> "$category"_significant_families_count.txt
grep -c "y" "$category"_singlelambda_k"$n"/*_family_results.txt >> "$category"_significant_families_count.txt
done
```

To write the just significant families to a file:
```
grep "y" Gamma_family_results.txt > Significant_families.txt
```

* proteases
```
category=proteasefamily
tree=/data/ccallen/2023_02_Sporothrix/09_orthology/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree
input=/data/ccallen/2023_02_Sporothrix/14_CAFE/protease_families_funannotate.txt


cd /data/ccallen/2023_02_Sporothrix/14_CAFE
for n in {1..10}
do
output="$category"_singlelambda_k"$n"
mkdir "$output"
/data/ccallen/bin/CAFE5/bin/cafe5 -t "$tree" -i "$input" -p -k "$n" -o "$output" &>> "$output"/cafe.log
done

for n in {1..10}
do
echo "$category"_singlelambda_k"$n" >> "$category"_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"/Base_results.txt >> "$category"_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"/Gamma_results.txt >> "$category"_discrete_rate_categories.txt
done

# count significant families at the p=0.05 threshold:
for n in {1..10}
do
echo "$category"_singlelambda_k"$n" >> "$category"_significant_families_count.txt
grep -c "y" "$category"_singlelambda_k"$n"/*_family_results.txt >> "$category"_significant_families_count.txt
done
```


* BGC (failed!)
```
category=secmets
tree=/data/ccallen/2023_02_Sporothrix/09_orthology/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree
input=/data/ccallen/2023_02_Sporothrix/14_CAFE/secmets.txt

cd /data/ccallen/2023_02_Sporothrix/14_CAFE
for n in {1..10}
do
output="$category"_singlelambda_k"$n"_2
mkdir "$output"
/data/ccallen/bin/CAFE5/bin/cafe5 -t "$tree" -i "$input" -p -k "$n" -o "$output" &>> "$output"/cafe.log
done

for n in {1..10}
do
echo "$category"_singlelambda_k"$n"_2 >> "$category"_2_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"_2/Base_results.txt >> "$category"_2_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"_2/Gamma_results.txt >> "$category"_2_discrete_rate_categories.txt
done

# count significant families at the p=0.05 threshold:
for n in {1..10}
do
echo "$category"_singlelambda_k"$n"_2 >> "$category"_2_significant_families_count.txt
grep -c "y" "$category"_singlelambda_k"$n"_2/*_family_results.txt >> "$category"_2_significant_families_count.txt
done
```

* orthogroups
```
category=orthogroups
tree=/data/ccallen/2023_02_Sporothrix/09_orthology/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.coresporothrix.contree
input=/data/ccallen/2023_02_Sporothrix/14_CAFE/orthogroups.txt


cd /data/ccallen/2023_02_Sporothrix/14_CAFE
for n in {1..10}
do
output="$category"_singlelambda_k"$n"
mkdir "$output"
/data/ccallen/bin/CAFE5/bin/cafe5 -t "$tree" -i "$input" -p -k "$n" -o "$output" &>> "$output"/cafe.log
done

for n in {1..10}
do
echo "$category"_singlelambda_k"$n" >> "$category"_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"/Base_results.txt >> "$category"_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"/Gamma_results.txt >> "$category"_discrete_rate_categories.txt
done

# count significant families at the p=0.05 threshold:
for n in {1..10}
do
echo "$category"_singlelambda_k"$n" >> "$category"_significant_families_count.txt
grep -c "y" "$category"_singlelambda_k"$n"/*_family_results.txt >> "$category"_significant_families_count.txt
done
```

* use CafePlotter v0.20 to plot CAFE5 gene family expansion/contraction results
```
cafeplotter -i cazymefamily_singlelambda_k1 -o cazymefamily_singlelambda_k1_cafeplotter --format png  --fig_width 20 --fig_height 0.4 --count_label_size 10
cafeplotter -i proteasefamily_singlelambda_k1 -o proteasefamily_singlelambda_k1_cafeplotter --format png  --fig_width 20 --fig_height 0.4 --count_label_size 10
cafeplotter -i secmets_singlelambda_k1 -o secmets_singlelambda_k1_cafeplotter --format png  --fig_width 20 --fig_height 0.4 --count_label_size 10
cafeplotter -i orthogroups_singlelambda_k1 -o orthogroups_singlelambda_k1_cafeplotter --format png  --fig_width 20 --fig_height 0.4 --count_label_size 10
```





