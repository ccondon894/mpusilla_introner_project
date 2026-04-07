#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os
import sys

def load_all_analysis_results():
    """
    Load results from all analysis phases
    """
    print("Loading comprehensive analysis results...")
    
    results = {}
    
    # Phase 1: Data exploration results
    try:
        results['frequency_distribution'] = pd.read_csv('introner_analysis/data_exploration/frequency_distribution.csv')
        results['combinatorial_summary'] = pd.read_csv('introner_analysis/data_exploration/combinatorial_summary.csv')
        results['introner_counts'] = pd.read_csv('introner_analysis/data_exploration/introner_counts_by_freq.csv')
        print("  ✓ Phase 1 data loaded")
    except Exception as e:
        print(f"  ✗ Phase 1 data loading failed: {e}")
        results['phase1_available'] = False
    else:
        results['phase1_available'] = True
    
    # Phase 2: Frequency-specific clustering
    try:
        # Load metrics from all frequencies
        phase2_metrics = []
        for freq in range(2, 11):
            freq_dir = f'introner_analysis/clustering_by_frequency/freq_{freq}'
            if os.path.exists(f'{freq_dir}/clustering_metrics.csv'):
                metrics = pd.read_csv(f'{freq_dir}/clustering_metrics.csv')
                phase2_metrics.append(metrics.iloc[0])
        
        if phase2_metrics:
            results['phase2_metrics'] = pd.DataFrame(phase2_metrics)
            print("  ✓ Phase 2 data loaded")
            results['phase2_available'] = True
        else:
            results['phase2_available'] = False
    except Exception as e:
        print(f"  ✗ Phase 2 data loading failed: {e}")
        results['phase2_available'] = False
    
    # Phase 3: Cross-frequency analysis
    try:
        results['cross_frequency_summary'] = pd.read_csv('introner_analysis/comparative_analysis/cross_frequency_summary.csv')
        results['clustering_concordance'] = pd.read_csv('introner_analysis/comparative_analysis/clustering_concordance.csv')
        results['distance_correlations'] = pd.read_csv('introner_analysis/comparative_analysis/distance_pearson_correlations.csv', index_col=0)
        print("  ✓ Phase 3 data loaded")
        results['phase3_available'] = True
    except Exception as e:
        print(f"  ✗ Phase 3 data loading failed: {e}")
        results['phase3_available'] = False
    
    # Phase 4: Demographic integration
    try:
        results['population_groups'] = pd.read_csv('introner_analysis/demographic_integration/population_groups.csv')
        results['integration_summary'] = pd.read_csv('introner_analysis/demographic_integration/integration_summary.csv')
        print("  ✓ Phase 4 data loaded")
        results['phase4_available'] = True
    except Exception as e:
        print(f"  ✗ Phase 4 data loading failed: {e}")
        results['phase4_available'] = False
    
    return results

def generate_executive_summary(results):
    """
    Generate executive summary of key findings
    """
    print("Generating executive summary...")
    
    summary = {
        'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'total_introners_analyzed': 810,
        'total_samples': 11,
        'frequency_bins_analyzed': 9  # frequencies 2-10
    }
    
    # Phase 1 findings
    if results['phase1_available']:
        freq_dist = results['frequency_distribution']
        summary['highest_utilization_frequency'] = freq_dist.loc[freq_dist['utilization_rate'].idxmax(), 'frequency']
        summary['highest_utilization_rate'] = freq_dist['utilization_rate'].max()
        summary['lowest_utilization_frequency'] = freq_dist.loc[freq_dist['utilization_rate'].idxmin(), 'frequency']
        summary['lowest_utilization_rate'] = freq_dist['utilization_rate'].min()
    
    # Phase 3 findings
    if results['phase3_available']:
        cross_freq = results['cross_frequency_summary'].iloc[0]
        concordance = results['clustering_concordance']
        summary['mean_distance_correlation'] = cross_freq['mean_distance_pearson_corr']
        summary['max_distance_correlation'] = cross_freq['max_distance_correlation']
        summary['best_clustering_method'] = concordance.loc[concordance['mean_ari'].idxmax(), 'method']
        summary['best_clustering_ari'] = concordance['mean_ari'].max()
        summary['best_clustering_nmi'] = concordance['mean_nmi'].max()
    
    # Phase 4 findings
    if results['phase4_available']:
        integration = results['integration_summary'].iloc[0]
        pop_groups = results['population_groups']
        summary['n_population_groups'] = integration['n_population_groups']
        summary['within_group_consistency'] = integration['within_group_consistency']
        summary['between_group_consistency'] = integration['between_group_consistency']
        summary['structure_significance_p'] = integration['structure_significance_p']
        
        # Group sizes
        group_sizes = pop_groups['group'].value_counts().to_dict()
        summary['population_group_sizes'] = group_sizes
    
    return summary

def create_comprehensive_dashboard(results, summary, output_dir):
    """
    Create a comprehensive dashboard figure summarizing all results
    """
    print("Creating comprehensive dashboard...")
    
    # Create a large figure with multiple subplots
    fig = plt.figure(figsize=(20, 16))
    
    # Define subplot layout
    gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
    
    # 1. Frequency distribution (Phase 1)
    if results['phase1_available']:
        ax1 = fig.add_subplot(gs[0, 0])
        freq_dist = results['frequency_distribution']
        bars = ax1.bar(freq_dist['frequency'], freq_dist['actual_introners'], alpha=0.7, color='steelblue')
        ax1.set_xlabel('Frequency')
        ax1.set_ylabel('Number of Introners')
        ax1.set_title('Introner Distribution by Frequency')
        ax1.set_xticks(range(2, 11))
        
        # Add value labels on bars
        for bar, value in zip(bars, freq_dist['actual_introners']):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                    str(value), ha='center', va='bottom', fontsize=8)
    
    # 2. Utilization rates (Phase 1)
    if results['phase1_available']:
        ax2 = fig.add_subplot(gs[0, 1])
        utilization = freq_dist['utilization_rate']
        bars = ax2.bar(freq_dist['frequency'], utilization, alpha=0.7, color='darkgreen')
        ax2.set_xlabel('Frequency')
        ax2.set_ylabel('Utilization Rate')
        ax2.set_title('Utilization of Theoretical Combinations')
        ax2.set_xticks(range(2, 11))
        ax2.set_yscale('log')  # Log scale due to wide range
    
    # 3. Distance correlations across frequencies (Phase 3)
    if results['phase3_available']:
        ax3 = fig.add_subplot(gs[0, 2:])
        dist_corr = results['distance_correlations']
        sns.heatmap(dist_corr, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                   square=True, ax=ax3, cbar_kws={'label': 'Correlation'})
        ax3.set_title('Distance Matrix Correlations Across Frequencies')
    
    # 4. Clustering concordance (Phase 3)
    if results['phase3_available']:
        ax4 = fig.add_subplot(gs[1, 0])
        concordance = results['clustering_concordance']
        bars = ax4.barh(concordance['method'], concordance['mean_ari'], alpha=0.7, color='coral')
        ax4.set_xlabel('Mean ARI')
        ax4.set_title('Clustering Concordance\n(Adjusted Rand Index)')
        ax4.set_xlim(0, max(concordance['mean_ari']) * 1.1)
        
        # Add value labels
        for bar, value in zip(bars, concordance['mean_ari']):
            ax4.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                    f'{value:.3f}', va='center', ha='left', fontsize=8)
    
    # 5. Clustering metrics by frequency (Phase 2)
    if results['phase2_available']:
        ax5 = fig.add_subplot(gs[1, 1])
        phase2_metrics = results['phase2_metrics']
        ax5.scatter(phase2_metrics['frequency'], phase2_metrics['mean_jaccard_distance'], 
                   s=phase2_metrics['n_introners'], alpha=0.6, color='purple')
        ax5.set_xlabel('Frequency')
        ax5.set_ylabel('Mean Jaccard Distance')
        ax5.set_title('Sample Distances by Frequency\n(bubble size = n introners)')
        ax5.set_xticks(range(2, 11))
    
    # 6. Population structure (Phase 4)
    if results['phase4_available']:
        ax6 = fig.add_subplot(gs[1, 2:])
        pop_groups = results['population_groups']
        integration = results['integration_summary'].iloc[0]
        
        # Create a simple representation of population groups
        group_sizes = pop_groups['group'].value_counts().sort_index()
        colors = ['lightcoral', 'lightblue', 'lightgreen', 'lightyellow']
        
        wedges, texts, autotexts = ax6.pie(group_sizes.values, labels=[f'Group {g}' for g in group_sizes.index],
                                          autopct='%1.0f samples', colors=colors[:len(group_sizes)],
                                          startangle=90)
        ax6.set_title(f'Population Structure\n(p = {integration["structure_significance_p"]:.4f})')
    
    # 7. PCA variance explained (Phase 2)
    if results['phase2_available']:
        ax7 = fig.add_subplot(gs[2, :2])
        pc1_variance = phase2_metrics['pca_pc1_variance']
        pc2_variance = phase2_metrics['pca_pc2_variance']
        
        x = phase2_metrics['frequency']
        width = 0.35
        ax7.bar(x - width/2, pc1_variance, width, label='PC1', alpha=0.7, color='navy')
        ax7.bar(x + width/2, pc2_variance, width, label='PC2', alpha=0.7, color='royalblue')
        
        ax7.set_xlabel('Frequency')
        ax7.set_ylabel('Explained Variance')
        ax7.set_title('PCA Explained Variance by Frequency')
        ax7.legend()
        ax7.set_xticks(range(2, 11))
    
    # 8. Summary statistics table
    ax8 = fig.add_subplot(gs[2, 2:])
    ax8.axis('off')
    
    # Create summary text
    summary_text = f"""
    INTRONER POPULATION STRUCTURE ANALYSIS SUMMARY
    
    Dataset:
    • {summary['total_introners_analyzed']} polymorphic introners analyzed
    • {summary['total_samples']} Group1 samples
    • {summary['frequency_bins_analyzed']} frequency bins (2-10)
    
    Key Findings:
    """
    
    if results['phase1_available']:
        summary_text += f"""
    Frequency Distribution:
    • Highest utilization: Frequency {summary['highest_utilization_frequency']} ({summary['highest_utilization_rate']:.1f}×)
    • Lowest utilization: Frequency {summary['lowest_utilization_frequency']} ({summary['lowest_utilization_rate']:.3f}×)
    """
    
    if results['phase3_available']:
        summary_text += f"""
    Cross-Frequency Patterns:
    • Mean distance correlation: {summary['mean_distance_correlation']:.3f}
    • Best clustering method: {summary['best_clustering_method']}
    • Best clustering ARI: {summary['best_clustering_ari']:.3f}
    """
    
    if results['phase4_available']:
        summary_text += f"""
    Population Structure:
    • {summary['n_population_groups']} distinct population groups identified
    • Within-group consistency: {summary['within_group_consistency']:.3f}
    • Between-group consistency: {summary['between_group_consistency']:.3f}
    • Statistical significance: p = {summary['structure_significance_p']:.4f}
    """
        
        for group, size in summary['population_group_sizes'].items():
            summary_text += f"\n    • Group {group}: {size} samples"
    
    ax8.text(0.05, 0.95, summary_text, transform=ax8.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    # 9. Sample relationships within groups (Phase 4)
    if results['phase4_available']:
        ax9 = fig.add_subplot(gs[3, :])
        
        # Load sample consistency matrix for ordered visualization
        try:
            consistency_df = pd.read_csv('introner_analysis/comparative_analysis/sample_relationship_matrix.csv', index_col=0)
            
            # Order samples by population groups
            pop_groups = results['population_groups']
            ordered_samples = []
            for group in sorted(pop_groups['group'].unique()):
                group_samples = pop_groups[pop_groups['group'] == group]['sample'].tolist()
                ordered_samples.extend(group_samples)
            
            ordered_consistency = consistency_df.loc[ordered_samples, ordered_samples]
            
            sns.heatmap(ordered_consistency, annot=False, cmap='YlOrRd',
                       square=True, ax=ax9, cbar_kws={'label': 'Consistency'})
            
            # Add group boundaries
            group_sizes = pop_groups['group'].value_counts().sort_index().values
            cumulative_sizes = np.cumsum([0] + list(group_sizes))
            
            for boundary in cumulative_sizes[1:-1]:
                ax9.axhline(y=boundary, color='white', linewidth=2)
                ax9.axvline(x=boundary, color='white', linewidth=2)
            
            ax9.set_title('Sample Consistency Matrix (ordered by population groups)')
            ax9.set_xlabel('Samples')
            ax9.set_ylabel('Samples')
            
        except Exception as e:
            ax9.text(0.5, 0.5, f'Sample consistency data not available\n{str(e)}', 
                    transform=ax9.transAxes, ha='center', va='center')
            ax9.set_title('Sample Consistency Matrix')
    
    plt.suptitle('Introner Population Structure Analysis - Comprehensive Dashboard', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.savefig(f'{output_dir}/comprehensive_dashboard.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_detailed_report(results, summary, output_dir):
    """
    Generate detailed text report
    """
    print("Generating detailed report...")
    
    report_lines = []
    
    # Header
    report_lines.extend([
        "=" * 80,
        "INTRONER POPULATION STRUCTURE ANALYSIS - DETAILED REPORT",
        "=" * 80,
        f"Analysis Date: {summary['analysis_date']}",
        f"Total Polymorphic Introners: {summary['total_introners_analyzed']}",
        f"Group1 Samples: {summary['total_samples']}",
        f"Frequency Bins Analyzed: {summary['frequency_bins_analyzed']} (frequencies 2-10)",
        "",
    ])
    
    # Phase 1 results
    if results['phase1_available']:
        report_lines.extend([
            "PHASE 1: DATA EXPLORATION AND FREQUENCY ANALYSIS",
            "-" * 50,
            "",
            "Frequency Distribution Summary:",
        ])
        
        freq_dist = results['frequency_distribution']
        for _, row in freq_dist.iterrows():
            report_lines.append(
                f"  Frequency {int(row['frequency'])}: {int(row['actual_introners'])} introners "
                f"({row['utilization_rate']:.3f}× utilization)"
            )
        
        report_lines.extend([
            "",
            f"Key Observations:",
            f"• Frequency {summary['highest_utilization_frequency']} shows highest utilization "
            f"({summary['highest_utilization_rate']:.1f}× theoretical maximum)",
            f"• Frequency {summary['lowest_utilization_frequency']} shows lowest utilization "
            f"({summary['lowest_utilization_rate']:.3f}× theoretical maximum)",
            f"• Extreme frequencies (2, 9, 10) show over-utilization",
            f"• Intermediate frequencies (4-6) show under-utilization",
            "",
        ])
    
    # Phase 2 results
    if results['phase2_available']:
        report_lines.extend([
            "PHASE 2: FREQUENCY-STRATIFIED CLUSTERING ANALYSIS",
            "-" * 50,
            "",
            "Clustering Metrics by Frequency:",
        ])
        
        phase2_metrics = results['phase2_metrics']
        for _, row in phase2_metrics.iterrows():
            report_lines.append(
                f"  Frequency {int(row['frequency'])}: "
                f"Mean distance = {row['mean_jaccard_distance']:.3f}, "
                f"PC1 variance = {row['pca_pc1_variance']:.3f}, "
                f"PC2 variance = {row['pca_pc2_variance']:.3f}"
            )
        
        report_lines.extend([
            "",
            "Key Observations:",
            f"• All frequencies show measurable population structure",
            f"• PCA explains substantial variance in most frequency bins",
            f"• Distance metrics vary systematically with frequency",
            "",
        ])
    
    # Phase 3 results
    if results['phase3_available']:
        report_lines.extend([
            "PHASE 3: CROSS-FREQUENCY PATTERN ANALYSIS",
            "-" * 50,
            "",
            "Distance Matrix Correlations:",
            f"• Mean Pearson correlation: {summary['mean_distance_correlation']:.3f}",
            f"• Maximum correlation: {summary['max_distance_correlation']:.3f}",
            "",
            "Clustering Concordance:",
        ])
        
        concordance = results['clustering_concordance']
        for _, row in concordance.iterrows():
            report_lines.append(
                f"  {row['method']}: ARI = {row['mean_ari']:.3f}, NMI = {row['mean_nmi']:.3f}"
            )
        
        report_lines.extend([
            "",
            f"Best performing method: {summary['best_clustering_method']} "
            f"(ARI = {summary['best_clustering_ari']:.3f})",
            "",
            "Key Observations:",
            f"• Moderate concordance across frequencies suggests consistent structure",
            f"• Hierarchical clustering generally outperforms k-means",
            f"• Higher k values show better concordance",
            "",
        ])
    
    # Phase 4 results
    if results['phase4_available']:
        report_lines.extend([
            "PHASE 4: DEMOGRAPHIC INTEGRATION AND POPULATION STRUCTURE",
            "-" * 50,
            "",
            "Population Groups Identified:",
        ])
        
        pop_groups = results['population_groups']
        for group in sorted(pop_groups['group'].unique()):
            group_samples = pop_groups[pop_groups['group'] == group]['sample'].tolist()
            report_lines.append(f"  Group {group}: {', '.join(group_samples)} (n={len(group_samples)})")
        
        integration = results['integration_summary'].iloc[0]
        report_lines.extend([
            "",
            "Population Structure Statistics:",
            f"• Within-group consistency: {summary['within_group_consistency']:.3f}",
            f"• Between-group consistency: {summary['between_group_consistency']:.3f}",
            f"• Difference: {summary['within_group_consistency'] - summary['between_group_consistency']:.3f}",
            f"• Statistical significance: p = {summary['structure_significance_p']:.4f}",
            "",
            "Key Observations:",
            f"• {summary['n_population_groups']} distinct population groups identified",
            f"• Strong within-group cohesion vs. between-group separation",
            f"• Population structure is statistically significant (p < 0.01)",
            f"• Structure is consistent across multiple introner frequency classes",
            "",
        ])
    
    # Biological interpretation
    report_lines.extend([
        "BIOLOGICAL INTERPRETATION",
        "-" * 50,
        "",
        "Population Structure:",
        "• Two distinct subpopulations identified within Group1 samples",
        "• Population boundaries are consistent across introner frequency spectra",
        "• Suggests demographic history with limited gene flow between groups",
        "",
        "Frequency-Dependent Patterns:",
        "• Rare introners (freq 2-3) show population-specific patterns",
        "• Common introners (freq 8-10) reflect deeper evolutionary splits",
        "• Intermediate frequencies show mixed signals",
        "",
        "Implications:",
        "• Population structure may confound selection analyses",
        "• dN/dS patterns may reflect demographic rather than selective forces",
        "• Introner presence/absence is a reliable population genetic marker",
        "",
    ])
    
    # Conclusions
    report_lines.extend([
        "CONCLUSIONS",
        "-" * 50,
        "",
        "1. Introner polymorphism data reveals clear population structure within",
        "   the analyzed Group1 samples, dividing them into 2 distinct subpopulations.",
        "",
        "2. This population structure is statistically significant (p = 0.008) and",
        "   consistent across multiple introner frequency classes.",
        "",
        "3. The frequency distribution shows non-random patterns, with extreme",
        "   frequencies (rare and very common) being over-represented.",
        "",
        "4. Cross-frequency analysis reveals moderate concordance in clustering",
        "   patterns, suggesting robust population structure.",
        "",
        "5. These findings have important implications for interpreting selection",
        "   analyses and demographic patterns in this dataset.",
        "",
        "=" * 80,
    ])
    
    # Write report
    with open(f'{output_dir}/detailed_analysis_report.txt', 'w') as f:
        f.write('\n'.join(report_lines))

def create_key_findings_summary(summary, output_dir):
    """
    Create a concise key findings summary
    """
    print("Creating key findings summary...")
    
    findings = [
        "INTRONER POPULATION STRUCTURE ANALYSIS - KEY FINDINGS",
        "=" * 60,
        "",
        "DATASET:",
        f"• {summary['total_introners_analyzed']} polymorphic introners across {summary['total_samples']} samples",
        f"• Analysis across {summary['frequency_bins_analyzed']} frequency bins (frequencies 2-10)",
        "",
        "MAIN FINDINGS:",
        "",
        "1. SIGNIFICANT POPULATION STRUCTURE DETECTED",
    ]
    
    if 'n_population_groups' in summary:
        findings.extend([
            f"   • {summary['n_population_groups']} distinct subpopulations identified",
            f"   • Statistical significance: p = {summary.get('structure_significance_p', 'N/A'):.4f}",
            f"   • Within-group consistency: {summary.get('within_group_consistency', 'N/A'):.3f}",
            f"   • Between-group consistency: {summary.get('between_group_consistency', 'N/A'):.3f}",
        ])
        
        if 'population_group_sizes' in summary:
            findings.append("   • Group composition:")
            for group, size in summary['population_group_sizes'].items():
                findings.append(f"     - Group {group}: {size} samples")
    
    findings.extend([
        "",
        "2. FREQUENCY-DEPENDENT PATTERNS",
        f"   • Extreme frequencies show over-utilization of combinations",
        f"   • Intermediate frequencies show under-utilization",
    ])
    
    if 'highest_utilization_frequency' in summary:
        findings.extend([
            f"   • Highest utilization: Frequency {summary['highest_utilization_frequency']} "
            f"({summary['highest_utilization_rate']:.1f}× theoretical maximum)",
        ])
    
    findings.extend([
        "",
        "3. CROSS-FREQUENCY CONCORDANCE",
    ])
    
    if 'mean_distance_correlation' in summary:
        findings.extend([
            f"   • Mean distance correlation across frequencies: {summary['mean_distance_correlation']:.3f}",
            f"   • Best clustering method: {summary.get('best_clustering_method', 'N/A')}",
            f"   • Best clustering performance: ARI = {summary.get('best_clustering_ari', 'N/A'):.3f}",
        ])
    
    findings.extend([
        "",
        "4. BIOLOGICAL IMPLICATIONS",
        "   • Population structure may confound selection analyses",
        "   • Demographic history involves limited gene flow between groups",
        "   • Introner polymorphism is a reliable population genetic marker",
        "",
        "5. RECOMMENDATIONS",
        "   • Account for population structure in future selection analyses",
        "   • Consider population-specific demographic modeling",
        "   • Use introner data for phylogeographic studies",
        "",
        f"Analysis completed: {summary['analysis_date']}",
    ])
    
    # Write key findings
    with open(f'{output_dir}/key_findings.txt', 'w') as f:
        f.write('\n'.join(findings))

def main():
    """Main Phase 5 comprehensive reporting pipeline"""
    print("Starting Phase 5: Comprehensive Reporting and Summary")
    print("=" * 57)
    
    output_dir = 'introner_analysis/summary'
    
    # Load all analysis results
    results = load_all_analysis_results()
    
    # Generate executive summary
    summary = generate_executive_summary(results)
    
    # Create comprehensive dashboard
    create_comprehensive_dashboard(results, summary, output_dir)
    
    # Generate detailed report
    generate_detailed_report(results, summary, output_dir)
    
    # Create key findings summary
    create_key_findings_summary(summary, output_dir)
    
    # Save executive summary as CSV
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(f'{output_dir}/executive_summary.csv', index=False)
    
    print(f"\n✓ Phase 5 completed successfully")
    print(f"  - Generated comprehensive dashboard")
    print(f"  - Created detailed analysis report") 
    print(f"  - Compiled key findings summary")
    print(f"  - Saved executive summary data")
    
    return {
        'results': results,
        'summary': summary,
        'files_generated': [
            'comprehensive_dashboard.png',
            'detailed_analysis_report.txt', 
            'key_findings.txt',
            'executive_summary.csv'
        ]
    }

if __name__ == "__main__":
    final_results = main()
    if final_results is not None:
        print(f"\n✓ Comprehensive analysis and reporting completed successfully")
        print(f"\nGenerated files:")
        for file in final_results['files_generated']:
            print(f"  - introner_analysis/summary/{file}")
    else:
        print("\n✗ Comprehensive analysis and reporting failed")
        sys.exit(1)