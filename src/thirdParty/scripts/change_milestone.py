#!/usr/bin/env python 
"""
Change the milestone (use the latest created) on all the open issues part of gitlab.com/PATO/PATO-dev.

Author: jeremie.meurisse@gmail.com
Gitlab Account: pato.devel.2@gmail.com
Usage: python change_milestone.py
"""

import selenium, time, sys
import undetected_chromedriver as uc
from selenium.webdriver.common.by import By
from getpass import getpass

username=str(input("Username: "))
password=str(getpass())

class MyUDC(uc.Chrome):
    def __del__(self):
        try:
            self.service.process.kill()
        except:  # noqa
            pass
driver = MyUDC()

# Access gitlab sign in
issues_gitlab_page="https://gitlab.com/PATO/PATO-dev/-/issues/?sort=created_date&state=opened&first_page_size=100"
driver.get(issues_gitlab_page)
time.sleep(3)

elem_tag_name='input'
elems=driver.find_elements(By.TAG_NAME,elem_tag_name)
found=False
for elem_i in elems:
    if elem_i.get_attribute("id") == "user_login":
        elem_i.send_keys(username)
        found=True
        break
if not found:
    print("Error: "+elem_tag_name+" not found.")
    sys.exit()
time.sleep(2)

elem_tag_name='input'
elems=driver.find_elements(By.TAG_NAME,elem_tag_name)
found=False
for elem_i in elems:
    if elem_i.get_attribute("id") == "user_password":
        elem_i.send_keys(password)
        found=True
        break
if not found:
    print("Error: "+elem_tag_name+" not found.")
    sys.exit()
time.sleep(2)

elem_tag_name='button'
elems=driver.find_elements(By.TAG_NAME,elem_tag_name)
found=False
for elem_i in elems:
    if elem_i.get_attribute("data-testid") == "sign-in-button":
        elem_i.click()
        found=True
        break
if not found:
    print("Error: "+elem_tag_name+" not found.")
    sys.exit()
time.sleep(2)

# Add the first milestone to all the issues
elem_tag_name='a'
elems=driver.find_elements(By.TAG_NAME,elem_tag_name) # all the issues
hrefs=[]
for elem_i in elems:
    if elem_i.get_attribute("dir") == "auto":
        href=elem_i.get_attribute("href")
        hrefs.append(href)

for href_i in hrefs:
    driver.get(href_i)
    time.sleep(4)
    elem_tag_name='button'
    elems=driver.find_elements(By.TAG_NAME,elem_tag_name)
    found=False
    for elem_i in elems:
        if elem_i.get_attribute("id") == "milestone-edit":
            elem_i.click()
            found=True
            break
    if not found:
        print("Error: "+elem_tag_name+" milestone-edit not found for "+href_i)
    time.sleep(4)
    elem_tag_name='button'
    elems=driver.find_elements(By.TAG_NAME,elem_tag_name)
    found=False
    for elem_i in elems:
        if elem_i.get_attribute("data-testid") == "milestone-items":
            elem_i.click()
            found=True
            break
    if not found:
        print("Error: "+elem_tag_name+" milestone-items not found for "+href_i)
    time.sleep(2)
