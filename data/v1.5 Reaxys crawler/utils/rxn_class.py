import re
import pandas as pd

class RXN():
    # rxn的自身参数
    def __init__(self):
        self.index = str() # 反应索引

        self.reactants = list() # 反应物
        self.products = list() # 产物

        self.reagents = list() # 试剂
        self.cats = list() # 催化剂
        self.solvents = list() # 溶剂
        self.temp = str() # 温度
        self.time = str() # 时间

        self.rxn_yield = str() # 产率

        self.ref = str() # 参考文献
        self.rxn_id = str() # Reaxys ID

    # rxn信息获取
    def get_info(self, html_parser, i):
        """
        :param html_parser: html 解析器
        :param i: 该反应类型中的第i条条件
        """
        # 反应索引
        self.index = html_parser.xpath("//span[@class='rx-element-index']/text()")[0]

        rxn_index = html_parser.xpath("//div[@class='rx-reactions-table__conditions__steps']")[i]
        # 若不存在stages-row
        rxn_index = rxn_index.xpath('./div[@class="stages-row"]')
        if len(rxn_index) == 0:
            return
        else:
            rxn_index = rxn_index[0]
        # 试剂
        rxn_r = rxn_index.xpath('./span[@class="stage-reagents"]')
        if len(rxn_r) == 0:
            pass
        else:
            for rea_html in rxn_r:
                rea = str(rea_html.xpath("string()"))
                if len(re.findall(".+[^;\xa0]", rea)) == 0:
                    self.reagents.append(rea)
                else:
                    rea = re.findall(".+[^;\xa0]", rea)[0] # html文本中有奇怪字符
                    self.reagents.append(rea)
        # 催化剂
        rxn_c = rxn_index.xpath("./span[@class='stage-catalyst']")
        if len(rxn_c) == 0:
            pass
        else:
            for cat_html in rxn_c:
                cat = str(cat_html.xpath("string()"))
                if len(re.findall(".+[^;\xa0]", cat)) == 0:
                    self.cats.append(cat)
                else:
                    cat = re.findall(".+[^;\xa0]", cat)[0] # html文本中有奇怪字符
                    self.cats.append(cat)
        # 溶剂
        rxn_s = rxn_index.xpath("./span[@class='stage-solvents']")
        if len(rxn_s) == 0:
            pass
        else:
            for sol_html in rxn_s:
                sol = str(sol_html.xpath("string()"))
                if len(re.findall(".+[^;\xa0]", sol)) == 0:
                    self.solvents.append(sol)
                else:
                    sol = re.findall(".+[^;\xa0]", sol)[0] # html文本中有奇怪字符
                    self.solvents.append(sol)
        # 时间 & 温度
        rxn_cond = rxn_index.xpath("string(./span[@class='conditions'])")
        temp = re.findall("(?<=at\s)-?[0-9]+\.[0-9]*|-?[0-9]+(?=℃)", rxn_cond)
        time = re.findall("(?<=for\s)-?[0-9]+\.[0-9]*|-?[0-9]+(?=h)", rxn_cond)
        if len(temp) == 0:
            self.temp = "/"
        else:
            self.temp = temp[0]

        if len(time) == 0:
            self.time = "/"
        else:
            self.time = time[0]

        # 产率
        rxn_y = html_parser.xpath("//td[@class='rx-reactions-table__yield display-table-cell']")[i]
        if len(rxn_y) == 0:
            self.rxn_yield = "/"
        else:
            self.rxn_yield = rxn_y.xpath("string()")

        # 参考文献
        rxn_ref = html_parser.xpath("//div[@class='citation clear']")[i]
        if len(rxn_ref) == 0:
            self.ref = "/"
        else:
            self.ref = rxn_ref.xpath("string()")

    def show_info(self):
        print(self.index, self.reactants, self.products, self.reagents, self.cats, self.solvents, self.temp, self.time, self.rxn_yield, self.rxn_id, self.ref)

# 将pandas的dataframe转为rxn的list
def df_to_rxn_list(df):
    """
    :param df: pandas的dataframe
    :return: 所有反应的rxn_class，以list类型返回
    """

    rxn_list = list()  # 反应的list
    data_size = df.shape[0]  # 获取反应数量

    for num in range(data_size):
        rxn = RXN() # 实例化反应的类
        # 获取 index
        rxn.index = df.loc[num]["rxn_index"]
        # 获取 temperature /C
        rxn.temp = df.loc[num]["temperature /C"]
        # 获取 time /h
        rxn.time = df.loc[num]["time /h"]
        # 获取 Yield
        rxn.rxn_yield = df.loc[num]["Yield"]
        # 获得Reaction ID
        rxn.rxn_id = df.loc[num]["Reaction ID"]
        # 获得参考文献
        rxn.ref = df.loc[num]["Reference"]

        # 获取 reactant
        for col in df.columns:
            if "reactants" in col:
                if df.loc[num][col] != "/": # 去除 /
                    rxn.reactants.append(df.loc[num][col])
        # 获取 product
        for col in df.columns:
            if "products" in col:
                if df.loc[num][col] != "/":  # 去除 /
                    rxn.products.append(df.loc[num][col])
        # 获取 reagent
        for col in df.columns:
            if "reagents" in col:
                if df.loc[num][col] != "/":  # 去除 /
                    rxn.reagents.append(df.loc[num][col])
        # 获取 catalyst
        for col in df.columns:
            if "catalysts" in col:
                if df.loc[num][col] != "/":  # 去除 /
                    rxn.cats.append(df.loc[num][col])
        # 获取 solvents
        for col in df.columns:
            if "solvents" in col:
                if df.loc[num][col] != "/":  # 去除 /
                    rxn.solvents.append(df.loc[num][col])

        rxn_list.append(rxn) # 放入list

    return rxn_list

# 将rxn_list转为df
def rxn_list_to_df(rxn_list):
# 遍历rxn_list, 找到reactants, products...等参数的最大数量,以确定dataframe的规格
    # 设置最大个数
    num_reactants = 0
    num_products = 0
    num_reagents = 0
    num_cats = 0
    num_sol = 0

    for rxn in rxn_list:
        num_reactants = max(len(rxn.reactants), num_reactants)
        num_products = max(len(rxn.products), num_products)
        num_reagents = max(len(rxn.reagents), num_reagents)
        num_cats = max(len(rxn.cats), num_cats)
        num_sol = max(len(rxn.solvents), num_sol)

    # 设置dataframe的索引columns
    cols = list()
    # 反应索引
    cols.append("rxn_index")
    # 反应物
    for i in range(num_reactants):
        cols.append("reactants %d" % (i + 1))
    # 产物
    for i in range(num_products):
        cols.append("products %d" % (i + 1))
    # 试剂
    for i in range(num_reagents):
        cols.append("reagents %d" % (i + 1))
    # 催化剂
    for i in range(num_cats):
        cols.append("catalysts %d" % (i + 1))
    # 溶剂
    for i in range(num_sol):
        cols.append("solvents %d" % (i + 1))
    # 温度
    cols.append("temperature /C")
    # 时间
    cols.append("time /h")
    # 产率
    cols.append("Yield")
    # Reaction ID
    cols.append("Reaction ID")
    # Citation
    cols.append("Reference")

    # 将数据导入dataframe
    data = list()
    for rxn in rxn_list:
        meta_data = list()
        # 索引
        meta_data.append(rxn.index)
        # 反应物
        for reactant in rxn.reactants:
            meta_data.append(reactant)
        while len(rxn.reactants) < num_reactants:
            meta_data.append("/")
            rxn.reactants.append("/")
        # 产物
        for product in rxn.products:
            meta_data.append(product)
        while len(rxn.products) < num_products:
            meta_data.append("/")
            rxn.products.append("/")
        # 试剂
        for reagent in rxn.reagents:
            meta_data.append(reagent)
        while len(rxn.reagents) < num_reagents:
            meta_data.append("/")
            rxn.reagents.append("/")
        # 催化剂
        for cat in rxn.cats:
            meta_data.append(cat)
        while len(rxn.cats) < num_cats:
            meta_data.append("/")
            rxn.cats.append("/")
        # 溶剂
        for sol in rxn.solvents:
            meta_data.append(sol)
        while len(rxn.solvents) < num_sol:
            meta_data.append("/")
            rxn.solvents.append("/")
        # 温度
        meta_data.append(rxn.temp)
        # 时间
        meta_data.append(rxn.time)
        # 产率
        meta_data.append(rxn.rxn_yield)
        # Reaction ID
        meta_data.append(rxn.rxn_id)
        # Reference
        meta_data.append(rxn.ref)

        # 将元数据放入数据中
        data.append(meta_data)

    # 生成dataframe
    df = pd.DataFrame(data, columns=cols)
    return df

if __name__ == "__main__":
    data = pd.read_excel("../data/raw data.xlsx")
    rxn_list = df_to_rxn_list(data)
    print(rxn_list[0].rxn_yield)
    print(rxn_list_to_df(rxn_list))