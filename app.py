import streamlit as st
import scanpy as sc
import plotly.express as px
import pandas as pd


st.set_page_config(layout="wide")

def load_sc_data(file_path):
    adata = sc.read_h5ad(file_path)
    return adata

def creat_umap_plot(adata,select:str):
    umap_df = pd.DataFrame(adata.obsm['X_umap'],columns=['UMAP1','UMAP2'])
    umap_df.index = adata.obs.index
    umap_df[select] = adata.obs[select].astype(str)
    fig = px.scatter(
        umap_df,
        x='UMAP1',y='UMAP2',
        color= select,
        labels=select,
        title =''
    )
    fig.update_yaxes(
        showgrid=False,
        zeroline=False
    )
    fig.update_layout(width=400,height=400)
    return fig

def creat_umap_exp_plot(adata,select):
    umap_df = pd.DataFrame(adata.obsm['X_umap'],columns=['UMAP1','UMAP2'])
    umap_df.index = adata.obs.index
    figs = []
    for gene in select:
        umap_df[gene] = adata[:, gene].X.toarray()

        fig = px.scatter(
            umap_df,
            x='UMAP1', y='UMAP2',
            color=gene,
            labels={gene: gene},
            title=f'UMAP Plot for {gene}'
        )
        fig.update_yaxes(
            showgrid=False,
            zeroline=False
        )
        fig.update_layout(width=300,height=300)
        figs.append(fig)
    
    return figs

def creat_heatmap_exp_plot(adata,gene_expression,select):
    # gene_expression = ['Sox17','St18']
    # select = 'cell_type'
    fig = sc.pl.matrixplot(adata, gene_expression, groupby=select, use_raw=True,return_fig=True)
    df = pd.DataFrame(fig.values_df)
    fig = px.imshow(df)
    return fig

def creat_dotplot_exp_plot(adata,gene_expression,select):
    fig = sc.pl.dotplot(adata, gene_expression, select, use_raw=False,return_fig=True)
    fig.dot_color_df['Y'] = fig.dot_color_df.index
    dot_color_df = fig.dot_color_df.melt(id_vars=['Y'])
    dot_size_df = fig.dot_size_df.melt()
    # 为 Bubble Plot 准备数据
    bubble_data = pd.DataFrame({
        'X': dot_color_df['variable'],
        'Y': dot_color_df['Y'],
        'Color': dot_color_df['value'],
        'Size': dot_size_df['value'],
        # 'Labels': dot_color_df['group']  # 根据需要修改标签
    })

    # 创建 Bubble Plot
    bubble_fig = px.scatter(bubble_data, x='X', y='Y', size='Size',color= 'Color')
    # # 设置图表布局
    bubble_fig.update_layout(
        title='Bubble Plot',
        xaxis_title='Gene',
        yaxis_title='Group'
    )
    bubble_fig.update_yaxes(
            showgrid=False,
            zeroline=False
        )
    return bubble_fig

def creat_violin_exp_plot(adata,gene_expression,select):
    fig = sc.pl.stacked_violin(adata,gene_expression, groupby=select,return_fig=True)
    fig.DEFAULT_CATEGORY_HEIGHT = 0.03*len(adata.obs[select].drop_duplicates().to_list())
    fig.DEFAULT_CATEGORY_WIDTH = 0.15*len(gene_expression)
    return fig

def get_diff_gene(adata):
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    all_diff = []
    for c in groups:
        tmp = sc.get.rank_genes_groups_df(adata, group=c)
        # tmp = tmp[(tmp['logfoldchanges'].abs()>1) & (tmp['pvals_adj']<0.05) & (tmp['pct_nz_group'] >= 0.2)]
        tmp['cluster'] = c
        all_diff.append(tmp)
    all_diff = pd.concat(all_diff, axis=0)
    all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    first_col = all_diff.pop('cluster')
    all_diff.insert(0, 'cluster', first_col)
    return all_diff

def main():
    st.title("单细胞浏览器")
    uploaded_file = st.file_uploader("上传单细胞数据文件（.h5ad格式）", type=["h5ad"], key="sc_data")
    if uploaded_file is not None:
        adata = load_sc_data(uploaded_file)
        # 添加选择框来选择基因表达
        st.sidebar.title('侧边栏')  # 侧边栏标题
        select = st.sidebar.selectbox('选择展示类别', ["choose an option"] + adata.obs.columns.tolist())
        gene_expression = st.sidebar.multiselect('选择基因表达', adata.var_names)
        col1 ,col2 = st.columns((2,4))
        if select != "choose an option":
            umap_fig = creat_umap_plot(adata,select)
            col1.plotly_chart(umap_fig, use_container_width=True)
        if gene_expression:
            # umap_fig = creat_umap_exp_plot(adata,gene_expression)
            # for fig in umap_fig:
            #     col2.plotly_chart(fig, use_container_width=True)
            with col2:
                tab1 ,tab2 ,tab3,tab4 = st.tabs(['UMAP','Bubble plot','Heatmap','Violin Plot'])
                with tab1:
                    num_cols = 3  # 设定你想要展示的列数
                    num_genes = len(gene_expression)
                    num_rows = -(-num_genes // num_cols)  # 向上取整得到行数

                    gene_column = st.columns(num_cols)  # 根据列数创建对应的列

                    for i, gene in enumerate(gene_expression):
                        col_idx = i % num_cols  # 当前列的索引
                        row_idx = i // num_cols  # 当前行的索引
                        fig = creat_umap_exp_plot(adata, [gene])
                        gene_column[col_idx].plotly_chart(fig[0], use_container_width=True)
                with tab2:
                    fig = creat_dotplot_exp_plot(adata,gene_expression,select)
                    st.plotly_chart(fig, use_container_width=True)
                with tab3:
                    fig = creat_heatmap_exp_plot(adata,gene_expression,select)
                    st.plotly_chart(fig, use_container_width=True)
                with tab4:
                    fig = creat_violin_exp_plot(adata,gene_expression,select)
                    st.pyplot(fig,use_container_width=False)

        # 在此处添加使用 Scanpy 和其他库对数据进行处理和可视化的代码
        with st.expander('Marker gene'):
            diff_genes = get_diff_gene(adata)
            selected_columns = st.selectbox('选择列',['no choose'] + diff_genes.iloc[:, 0].drop_duplicates().to_list())
            if selected_columns != 'no choose':
                diff_genes_select = diff_genes[diff_genes.iloc[:, 0] == selected_columns]
                st.dataframe(diff_genes_select,hide_index=True)
            else:
                st.dataframe(diff_genes,hide_index=True)


if __name__ == "__main__":
    main()
