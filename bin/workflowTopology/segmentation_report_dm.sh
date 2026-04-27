#!/usr/bin/env bash
while getopts "n:o:i:" OPTION; do
    case $OPTION in
        n) name="$OPTARG" ;;
        o) output_dir=$OPTARG ;;
        i) plot_statistics=$OPTARG ;;
    esac
done

output_file="${output_dir}/${name}.html"
[[ "${plot_statistics}" != "False" && "${plot_statistics}" != "false" ]] && output_file=$(realpath -s "$output_file")

STYLE="<style>body{font-family:sans-serif;background:#f4f7f6;padding:20px}.container{max-width:1000px;margin:0 auto;background:white;padding:30px;border-radius:8px;box-shadow:0 2px 10px rgba(0,0,0,0.1)}h1,h2{text-align:center;color:#2c3e50}.section{margin-bottom:30px;text-align:center}.grid{display:flex;flex-wrap:wrap;justify-content:center;gap:15px}.card{flex:1;min-width:280px;background:#fafafa;padding:10px;border:1px solid #eee;border-radius:4px}.label{font-weight:bold;display:block;margin-bottom:5px;color:#666}img{max-width:100%;height:auto;border-radius:4px}hr{border:0;border-top:1px dashed #ccc;margin:30px 0}</style>"

cat << EOF > "$output_file"
<!DOCTYPE html><html><head><title>DM Report - $name</title>$STYLE</head>
<body><div class="container"><h1>EpiSegMix Report: $name</h1>
<div class="section"><h2>Segmentation</h2><span class="label">State Colors</span><img src="${name}-state-colors.png"></div>
<div class="section"><h2>Normalized Counts</h2><img src="${name}-normEmission.png"></div>
<div class="section"><h2>Characteristics</h2><div class="grid">
<div class="card"><span class="label">Transition Matrix</span><img src="${name}-transitionMatrix.png"></div>
<div class="card"><span class="label">Average Length</span><img src="${name}-stateLength.png"></div>
<div class="card"><span class="label">Coverage</span><img src="${name}-stateMembership.png"></div>
</div></div>
EOF

if [[ "${plot_statistics}" != "False" && "${plot_statistics}" != "false" ]]; then
cat << EOF >> "$output_file"
<hr><div class="section"><h2>Emission Distribution</h2><img src="${name}-stateDistribution.png"></div>
<div class="section"><h2>Input Characteristics</h2><div class="grid">
<div class="card"><span class="label">Correlation</span><img src="${name}-correlation.png"></div>
<div class="card"><span class="label">Distributions</span><img src="${name}-histogram.png"></div>
</div></div>
EOF
fi
echo "</div></body></html>" >> "$output_file"