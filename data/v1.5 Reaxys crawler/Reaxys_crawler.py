from selenium.webdriver.common.keys import Keys
from selenium import webdriver
from time import sleep
from utils.rxn_class import *
import pandas as pd
import pyperclip
from bs4 import BeautifulSoup
from lxml import etree
from queue import Queue
import random

# 1.模拟登录
path = "geckodriver.exe"
bro = webdriver.Firefox(executable_path=path) # 实例化一个浏览器对象

# 1.1 Reaxys->大工界面
# 打开网站
bro.get("https://www.reaxys.com/#/institution")
sleep(5)
# 接受cookie
accept_cookie_btn = bro.find_element_by_id("onetrust-accept-btn-handler")
accept_cookie_btn.click()
# 输入学校信息
institution_search = bro.find_element_by_id("id-institution-search-input")
institution = input("Your institute:")
institution_search.send_keys(institution)
sleep(3)

chk = False
while chk:
    if input("Please log in through your institution. After logging in, enter OK to continue executing the program").lower() == "ok":
        chk = True

# 2.从 Reaxys 爬取数据
# 2.1搜索反应
for i in range(5):
    if "#/search/quick" in bro.current_url:
        bro.get("%s" % bro.current_url)
        break
    sleep(5)
    if i == 5:
        raise RuntimeError("Failed to go to the View Result interface")
sleep(3)
# 输入反应名称
q = bro.find_element_by_css_selector('[aria-describedby="new-input-validation-message"]')
q_input = input("请输入要查询的内容")
q.send_keys("%s" % q_input)
sleep(3)
# 点击搜索
find_btn = bro.find_element_by_css_selector('[ng-class="{\'rx-btn--disabled\': !buildQuery.query.text.length && !buildQuery.structurePreview}"]')
find_btn.click()
sleep(3)
# 点击 Reaction 的 View Result 按钮
while True:
    if len(bro.find_elements_by_css_selector('[rx-focus-first-element=".e2e-view-results"][aria-label="View reactions results for "]')) != 0:
        rxn_btn = bro.find_element_by_css_selector('[rx-focus-first-element=".e2e-view-results"][aria-label="View reactions results for "]')
        rxn_btn.click()
        break

# 2.2爬取数据
for i in range(5):
    if "#/results/reactions" in bro.current_url:
        bro.get("%s" % bro.current_url)
        break
    sleep(5)
    if i == 5:
        raise RuntimeError("Failed to go to the View Result interface")
# 设置网页加载时间
t = int(input("Please evaluate the website loading time (s) based on your current internet speed and recommend 20 seconds"))

# rxn_list 用于存储各反应
rxn_list = list()

# 设置爬取的目标页数
target_page = -1
while target_page <= 0:
    print("Please manually go to the page you need and enter the next question requirements")
    ans = input("Please enter the number of pages to be crawled in the future. If you want to crawl all pages, please enter ‘all’")
    if ans == "all":
        target_page = 99999
    else:
        target_page = int(ans)
cnt = 0

# 获取网页信息，直到没有“下一页”按钮
sleep(t)

while True:

    # 获取smi及其他信息
    html = bro.page_source
    smi = Queue()
    view_options_btns = bro.find_elements_by_css_selector(
        "[type='button'][class='els-btn els-btn-link els-btn-sm structure-options-button e2e-dropdown e2e-options ']")
    copy_btns = bro.find_elements_by_css_selector(
        "[type='button'][class='els-btn els-btn-link els-btn-sm e2e-drop-item e2e-copy-smiles']")
    num_vo_btns = len(view_options_btns)

    # 若view按钮与copy按钮数量不符，则说明部分化合物无法获取smi，搜索下一页
    if num_vo_btns != len(copy_btns):
        print("There are compounds on this page that cannot obtain SMILES")
        # 若没有下一页的按钮
        if len(bro.find_elements_by_css_selector('[aria-label="Next page"]')) == 0:
            break
        # 翻到下一页,并让网页加载t秒
        else:
            next_page_btn = bro.find_elements_by_css_selector('[aria-label="Next page"]')
            next_page_btn[0].send_keys(Keys.ENTER)
            sleep(t)
            continue

    # 依次按按钮，获取smi
    for i in range(num_vo_btns):
        view_options_btns[i].send_keys(Keys.ENTER)
        sleep(4)
        # 按钮被隐藏
        bro.execute_script("arguments[0].click();", copy_btns[i])
        sleep(4)
        smi_text = pyperclip.paste()
        smi.put(smi_text)

    # 若因网络问题导致无法获取所有化合物smi
    if smi.qsize() < num_vo_btns:
        print("Due to slow internet speed, it is not possible to retrieve all Copy SMILES bottom information for this page")
        print("This page is the %d page from now on" % cnt)

    # 输出时间
    print("Processing: %.2f" % (cnt / target_page * 100), "%")
    print("The current time required to crawl a page is %.2f min" % ((t + smi.qsize() * 6) / 60))
    print("Still needed: %.2f h" % ((t + smi.qsize() * 6) * (target_page - cnt) / 3600))

    # 解析原始 html
    soup = BeautifulSoup(html, "html.parser")
    rxns_html = soup.find_all("li", class_="rx-list-item")

    # 获取各类反应 html
    for rxn_html in rxns_html:
        # 定义html解析器
        html_parser = etree.HTML(str(rxn_html))
        # 确定rxn_id
        rxn_id = html_parser.xpath("//span[@class='rx-reaction-id']/text()")
        rxn_id = re.findall("\d+", rxn_id[0])[0]
        # 确定反应物，产物个数
        num_react = len(html_parser.xpath(
            "//div[@class='substances-item-wrapper substances-item-wrapper--reactant substance-image-small substance-image-container']"))
        num_prod = len(html_parser.xpath(
            "//div[@class='substances-item-wrapper substances-item-wrapper--product substance-image-small substance-image-container']"))
        # 获取各类反应反应物与产物的smi
        rxn_react = list()
        rxn_prod = list()
        for i in range(num_react):
            rxn_react.append(smi.get())
        for i in range(num_prod):
            rxn_prod.append(smi.get())
        # 获取各类反应的反应条件
        num = len(html_parser.xpath("//div[@class='rx-reactions-table__conditions__steps']"))
        for i in range(num):
            rxn = RXN()
            rxn.reactants = rxn_react
            rxn.products = rxn_prod
            rxn.rxn_id = rxn_id
            # 获取反应
            rxn.get_info(html_parser, i)
            rxn_list.append(rxn)

    # 若因网络问题导致有smi没有被记录
    if smi.empty() == False:
        print("Due to the bad internet, it is not possible to access all SMIELS on this page")
        print("This page is the %d page from now on" % cnt)

    # 对爬取的网页计数
    cnt += 1
    if cnt == target_page:
        print("Crawler ends normally")
        break

    if cnt % 5 == 0: # 每5页记录一次
        df = rxn_list_to_df(rxn_list)
        df.to_excel("./data/%s raw data.xlsx" % q_input)

    # 判断是否为尾页
    if len(bro.find_elements_by_css_selector('[aria-label="Next page"]')) == 0:
        print("It has reached the end page, and the crawler has ended normally")
        break
    else:
        # 翻到下一页,并让网页加载t秒
        next_page_btn = bro.find_elements_by_css_selector('[aria-label="Next page"]')
        next_page_btn[0].send_keys(Keys.ENTER)
        sleep(t + random.uniform(-5, 5))

sleep(5)
bro.quit()

df = rxn_list_to_df(rxn_list)
df.to_excel("./data/%s raw data.xlsx" % q_input)