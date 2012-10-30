# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant::Config.run do |config|
  config.vm.box = "lucid64"
  config.vm.box_url = "http://files.vagrantup.com/lucid64.box"
  config.vm.customize ["modifyvm", :id, "--memory", 2048]
  config.vm.share_folder "gemtools", "gemtools", "."
  config.vm.share_folder "gemtools-examples", "gemtools-examples", "../gemtools-examples"
  # Enable and configure the chef solo provisioner
  config.vm.provision :chef_solo do |chef|
    # We're going to download our cookbooks from the web
    chef.cookbooks_path = "chef/cookbooks"
    chef.add_recipe("basebox")
  end
end
