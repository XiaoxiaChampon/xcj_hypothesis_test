#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 07:09:57 2023

@author: Xiaoxia Champon and Chathura Jayalath

Purpose: Find users that tweets more than 10 times, tweets was retweeted,liked a lot from brandwatch data (Only Twitter source)
"""

#%%load the packages
import pandas as pd

import datetime

import glob 

import os.path

print(pd.__version__)
#%%functions and variables needed
twitter_measurements = ["Twitter Retweets", "Twitter Reply Count", "Twitter Likes"]
# gropuby based on the above
def get_measures(gdf):
    author_measurements = {}
    author_measurements["num_tweets"] = gdf[twitter_measurements].shape[0]
    author_measurements.update( gdf[twitter_measurements].sum().to_dict() )
    author_measurements["retweets_per_tweet"] = author_measurements["Twitter Retweets"] / author_measurements["num_tweets"]
    author_measurements["replies_per_tweet"] = author_measurements["Twitter Reply Count"] / author_measurements["num_tweets"]
    author_measurements["likes_per_tweet"] = author_measurements["Twitter Likes"] / author_measurements["num_tweets"]
    return pd.Series(author_measurements)
def get_author_union(path_dir, event_date, num_before,num_after, num_top_user, file_name_key="*"):
    files_list = glob.glob(os.path.join(path_dir, "*{}*.csv.zip".format(file_name_key)))
    event_df = pd.concat( [ pd.read_csv(f,skiprows=5,parse_dates=["Date"],low_memory=False) for f in files_list ] )
    start_date = event_date - datetime.timedelta(days = num_before)
    end_date = event_date + datetime.timedelta(days = num_after)
    index_to_drop = event_df[~(~event_df.Author.isna() & (start_date <= event_df["Date"]) & (event_df["Date"] <= end_date) & (event_df.Domain=="twitter.com"))].index
    event_df.drop(index = index_to_drop, inplace=True)
    author_data = event_df.groupby("Author").apply(lambda gdf: get_measures(gdf))
    top_authors = {measure : author_data[measure].sort_values(ascending=False)[:num_top_user] 
                   for measure in ["retweets_per_tweet","replies_per_tweet","likes_per_tweet","num_tweets"]}
    union_top_authors = list(set.union( *[set(v.index) for v in top_authors.values()] ))
    
    #author_list = list(event_df_new.Author.value_counts()[:num_top_user].index)
    author_data = author_data[author_data.index.isin(union_top_authors)]
    event_df = event_df[event_df.Author.isin(union_top_authors)][["Date","Author","Sentiment"]]
    return (event_df, author_data)
#%%
event_date = datetime.datetime(2022, 3, 16)
#%% Nov 11, 2023 get union of by reply, retweet, likes, tweets  #260 users
isarel_event_data, isarel_author_data= get_author_union('/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/data/isarel/', event_date, 1, 20, 100)
#%%
isarel_event_data
#%%
isarel_author_data



